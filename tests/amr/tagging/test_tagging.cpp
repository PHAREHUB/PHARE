
#include <cmath>
#include <algorithm>

#include "phare/phare.hpp"
#include "amr/tagging/tagger.hpp"
#include "amr/tagging/tagger_factory.hpp"
#include "amr/resources_manager/resources_manager.hpp"

#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/utilities/span.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"

using namespace PHARE::amr;



TEST(test_tagger, fromFactory)
{
    using phare_types = PHARE::PHARE_Types<1, 1, 2>;
    auto hybridTagger = TaggerFactory<phare_types>::make("HybridModel", "default");
    EXPECT_TRUE(hybridTagger != nullptr);

    auto badTagger = TaggerFactory<phare_types>::make("invalidModel", "invalidStrat");
    EXPECT_TRUE(badTagger == nullptr);
}



using Param   = std::vector<double>;
using RetType = std::shared_ptr<PHARE::core::Span<double>>;

RetType step1(Param const& x)
{
    std::vector<double> values(x.size());
    std::transform(std::begin(x), std::end(x), std::begin(values),
                   [](auto xx) { return std::tanh((xx - 0.52) / 0.05); });
    return std::make_shared<PHARE::core::VectorSpan<double>>(std::move(values));
}

RetType step2(Param const& x, Param const& y)
{
    throw std::runtime_error("fix me");
}

template<std::size_t dim>
auto constexpr step_fn()
{
    if constexpr (dim == 1)
        return &step1;
    if constexpr (dim == 2)
        return &step2;
}


template<std::size_t dim>
PHARE::initializer::PHAREDict createDict()
{
    using InitFunctionT        = PHARE::initializer::InitFunction<dim>;
    auto static constexpr step = step_fn<dim>();

    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"] = int{2};
    dict["ions"]["pop0"]["name"]   = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]   = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(step);


    dict["ions"]["pop0"]["particle_initializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop0"]["particle_initializer"]["charge"]         = -1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

    dict["ions"]["pop1"]["name"] = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"] = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(step);


    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(step);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(step);


    dict["ions"]["pop1"]["particle_initializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop1"]["particle_initializer"]["charge"]         = -1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

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



template<std::size_t dim_, std::size_t interp_, std::size_t refinedPartNbr_>
struct TaggingTestInfo
{
    auto static constexpr dim            = dim_;
    auto static constexpr interp         = interp_;
    auto static constexpr refinedPartNbr = refinedPartNbr_;
};


template<typename TaggingTestInfo_t>
struct TestTagger : public ::testing::Test
{
    auto static constexpr dim            = TaggingTestInfo_t::dim;
    auto static constexpr interp_order   = TaggingTestInfo_t::interp;
    auto static constexpr refinedPartNbr = TaggingTestInfo_t::refinedPartNbr;

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
        : layout{TestGridLayout<GridLayoutT>::make(20)}
        , Bx{"Bx", HybridQuantity::Scalar::Bx, layout.allocSize(HybridQuantity::Scalar::Bx)}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Bz{"Bz", HybridQuantity::Scalar::Bz, layout.allocSize(HybridQuantity::Scalar::Bz)}
        , Ex{"Ex", HybridQuantity::Scalar::Ex, layout.allocSize(HybridQuantity::Scalar::Ex)}
        , Ey{"Ey", HybridQuantity::Scalar::Ey, layout.allocSize(HybridQuantity::Scalar::Ey)}
        , Ez{"Ez", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
        , model{createDict<dim>()}
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

using TaggingTestInfos = testing::Types<TaggingTestInfo<1, 1, 2> /*, TaggingTestInfo<2, 1, 4>*/>;
TYPED_TEST_SUITE(TestTagger, TaggingTestInfos);

// TODOmaybe find a way to test the tagging?
TYPED_TEST(TestTagger, scaledAvg)
{
    /*
      auto strat = DefaultHybridTaggerStrategy<SinglePatchHybridModel>();
      strat.tag(model, layout, tags.data());
      {
          auto start
              = layout.physicalStartIndex(PHARE::core::QtyCentering::dual,
      PHARE::core::Direction::X); auto end =
      layout.physicalEndIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

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
  */
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int testResult = RUN_ALL_TESTS();
    return testResult;
}
