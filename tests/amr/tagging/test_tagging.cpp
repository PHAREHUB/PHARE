

#include "phare/phare.h"
#include "amr/tagging/tagger.h"
#include "amr/tagging/tagger_factory.h"
#include "amr/resources_manager/resources_manager.h"

#include "core/data/ndarray/ndarray_vector.h"
#include "core/models/hybrid_state.h"
#include "core/utilities/span.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>

using namespace PHARE::amr;




TEST(test_tagger, fromFactory)
{
    using phare_types = PHARE::PHARE_Types<1, 1, 2>;
    auto hybridTagger = TaggerFactory<phare_types>::make("HybridModel", "avg");
    EXPECT_TRUE(hybridTagger != nullptr);

    auto badTagger = TaggerFactory<phare_types>::make("invalidModel", "scaledAvg");
    EXPECT_TRUE(badTagger == nullptr);
}



using Param   = std::vector<double> const;
using RetType = std::shared_ptr<PHARE::core::Span<double>>;

RetType step(Param& x)
{
    std::vector<double> values(x.size());
    std::transform(std::begin(x), std::end(x), std::begin(values),
                   [](auto xx) { return std::tanh((xx - 0.52) / 0.05); });
    return std::make_shared<PHARE::core::VectorSpan<double>>(std::move(values));
}


using InitFunctionT = PHARE::initializer::InitFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"] = int{2};
    dict["ions"]["pop0"]["name"]   = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]   = 1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop0"]["ParticleInitializer"]["density"] = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(step);


    dict["ions"]["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop0"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

    dict["ions"]["pop1"]["name"] = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"] = 1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop1"]["ParticleInitializer"]["density"] = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(step);


    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(step);


    dict["ions"]["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop1"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(step);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(step);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(step);

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = 0.12;

    return dict;
}




struct TestTagger : public ::testing::Test
{
    static constexpr auto dim            = 1;
    static constexpr auto interp_order   = 1;
    static constexpr auto refinedPartNbr = 2;

    using phare_types = PHARE::PHARE_Types<dim, interp_order, refinedPartNbr>;
    using Electromag  = typename phare_types::Electromag_t;
    using Ions        = typename phare_types::Ions_t;
    using Electrons   = typename phare_types::Electrons_t;
    using Field       = typename phare_types::Field_t;
    using GridLayoutT = GridLayout<GridLayoutImplYee<dim, interp_order>>;

    struct SinglePatchHybridModel
    {
        using gridlayout_type           = GridLayout<GridLayoutImplYee<dim, interp_order>>;
        static auto constexpr dimension = dim;
        HybridState<Electromag, Ions, Electrons> state;
    };


    GridLayoutT layout;
    Field Bx, By, Bz;
    Field Ex, Ey, Ez;
    SinglePatchHybridModel model;
    std::vector<int> tags;

    TestTagger()
        : layout{{0.05}, {20}, {0.}}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Ey{"Ey", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Ez{"Ez", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , model{createDict()}
        , tags(20 + layout.nbrGhosts(PHARE::core::QtyCentering::dual))
    {
        model.state.electromag.E.setBuffer("EM_E_x", &Ex);
        model.state.electromag.E.setBuffer("EM_E_y", &Ey);
        model.state.electromag.E.setBuffer("EM_E_z", &Ez);
        model.state.electromag.B.setBuffer("EM_B_x", &Bx);
        model.state.electromag.B.setBuffer("EM_B_y", &By);
        model.state.electromag.B.setBuffer("EM_B_z", &Bz);
        model.state.electromag.initialize(layout);
    }
};


TEST_F(TestTagger, scaledAvg)
{
    auto strat = ScaledAvgHybridTaggerStrategy<SinglePatchHybridModel>();
    strat.tag(model, layout, tags.data());
    {
        auto start
            = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);
        auto end
            = layout.physicalEndIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

        auto endCell     = layout.nbrCells()[0] - 1;
        double threshold = 0.1;


        for (auto iCell = 0u, ix = start; iCell <= endCell; ++ix, ++iCell)
        {
            auto Bxavg = (Bx(ix - 1) + Bx(ix) + Bx(ix + 1)) / 3.;
            auto Byavg = (By(ix - 1) + By(ix) + By(ix + 1)) / 3.;
            auto Bzavg = (Bz(ix - 1) + Bz(ix) + Bz(ix + 1)) / 3.;

            auto diffx = std::abs(Bxavg - Bx(ix));
            auto diffy = std::abs(Byavg - By(ix));
            auto diffz = std::abs(Bzavg - Bz(ix));

            auto max = std::max({diffx, diffy, diffz});
            if (max > threshold)
            {
                tags[iCell] = 1;
            }
            else
                tags[iCell] = 0;
        }
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int testResult = RUN_ALL_TESTS();
    return testResult;
}
