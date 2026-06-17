


#include "simulator/simulator.hpp"
#include "amr/tagging/concrete_tagger.hpp"
#include "amr/tagging/tagging_criteria.hpp"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>
#include <algorithm>

using namespace PHARE::amr;

// dict in the runtime contract: method + nbr_quantities + Q{i}/{name,threshold}
PHARE::initializer::PHAREDict taggingDict(std::string const& method,
                                          std::vector<std::pair<std::string, double>> const& qtys)
{
    PHARE::initializer::PHAREDict dict;
    dict["method"]         = method;
    dict["nbr_quantities"] = static_cast<int>(qtys.size());
    for (std::size_t i = 0; i < qtys.size(); ++i)
    {
        auto const path        = "Q" + std::to_string(i);
        dict[path]["name"]      = qtys[i].first;
        dict[path]["threshold"] = qtys[i].second;
    }
    return dict;
}

TEST(test_tagger, constructsWithValidMethodAndQuantity)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::HybridModel_t;
    auto dict                  = taggingDict("default", {{"B", 0.2}});
    EXPECT_NO_THROW((ConcreteTagger<hybrid_model>{dict}));
    auto dictLohner = taggingDict("lohner", {{"B", 0.2}});
    EXPECT_NO_THROW((ConcreteTagger<hybrid_model>{dictLohner}));
}

TEST(test_tagger, throwsOnUnknownMethod)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::HybridModel_t;
    auto dict                  = taggingDict("invalidStrat", {{"B", 0.2}});
    EXPECT_THROW((ConcreteTagger<hybrid_model>{dict}), std::runtime_error);
}

TEST(test_tagger, throwsOnHybridRho)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::HybridModel_t;
    auto dict                  = taggingDict("default", {{"rho", 0.1}});
    EXPECT_THROW((ConcreteTagger<hybrid_model>{dict}), std::runtime_error);
}


// minimal 1D field for testing the pure criteria functions directly
struct MockField1D
{
    std::vector<double> data;
    double operator()(std::uint32_t i) const { return data[i]; }
};

TEST(test_criteria, lohnerSeparatesStepFromFlat)
{
    // tanh step centered at cell 20 over 40 cells
    std::size_t const n = 40;
    MockField1D f;
    f.data.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        f.data[i] = std::tanh((static_cast<double>(i) - 20.) / 1.5);

    std::vector<MockField1D const*> comps{&f};

    // Sample a flank of the step (cell 18), NOT the center: lohner is a second-
    // difference estimator and vanishes at the inflection point of a symmetric
    // tanh (a_p = -a_m, a_0 = 0). Curvature, and thus the indicator, is maximal
    // on the flanks.
    auto const atStep = lohnerIndicator<1>(comps, {18u});
    auto const atFlat = lohnerIndicator<1>(comps, {5u});
    EXPECT_GT(atStep, atFlat);
    EXPECT_GT(atStep, 10. * atFlat);
}

TEST(test_criteria, defaultSeparatesStepFromFlat)
{
    std::size_t const n = 40;
    MockField1D f;
    f.data.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        f.data[i] = std::tanh((static_cast<double>(i) - 20.) / 1.5);

    std::vector<MockField1D const*> comps{&f};

    auto const atStep = defaultIndicator<1>(comps, {20u});
    auto const atFlat = defaultIndicator<1>(comps, {5u});
    EXPECT_GT(atStep, atFlat);
}

TEST(test_criteria, parseTaggingMethod)
{
    EXPECT_EQ(parseTaggingMethod("default"), TaggingMethod::Default);
    EXPECT_EQ(parseTaggingMethod("lohner"), TaggingMethod::Lohner);
    EXPECT_THROW(parseTaggingMethod("nope"), std::runtime_error);
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
    dict["ions"]["nbrPopulations"] = std::size_t{2};
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
    auto static constexpr opts           = PHARE::SimOpts{dim, interp_order, refinedPartNbr};

    using phare_types = PHARE::solver::PHARE_Types<opts>;
    using Electromag  = phare_types::Electromag_t;
    using Ions        = phare_types::Ions_t;
    using Electrons   = phare_types::Electrons_t;
    using GridLayoutT = GridLayout<GridLayoutImplYee<dim, interp_order>>;

    struct SinglePatchHybridModel
    {
        using gridlayout_type           = GridLayout<GridLayoutImplYee<dim, interp_order>>;
        static auto constexpr dimension = dim;
        HybridState<Electromag, Ions, Electrons> state;
    };

    GridLayoutT layout;

    UsableVecField<dim> B, E;

    SinglePatchHybridModel model;
    std::vector<int> tags;

    TestTagger()
        : layout{TestGridLayout<GridLayoutT>::make(20)}
        , B{"EM_B", layout, HybridQuantity::Vector::B}
        , E{"EM_E", layout, HybridQuantity::Vector::E}
        , model{createDict<dim>()}
        , tags(20 + layout.nbrGhosts(PHARE::core::QtyCentering::dual))
    {
        B.set_on(model.state.electromag.B);
        E.set_on(model.state.electromag.E);
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
