

#include "phare_solver.hpp"
#include "simulator/simulator_def.hpp"

#include "amr/tagging/tagger_factory.hpp"
#include "amr/tagging/tagging_quantities.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

using namespace PHARE::amr;

TEST(test_tagger, fromFactoryValid)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::Hybrid::Model_t;
    PHARE::initializer::PHAREDict dict;
    dict["hybrid_method"] = std::string{"default"};
    dict["threshold"]     = 0.2;
    auto hybridTagger     = TaggerFactory<hybrid_model>::make(dict);
    EXPECT_TRUE(hybridTagger != nullptr);
}

TEST(test_tagger, fromFactoryInvalid)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::Hybrid::Model_t;
    PHARE::initializer::PHAREDict dict;
    dict["hybrid_method"] = std::string{"invalidStrat"};
    auto hybridTagger     = TaggerFactory<hybrid_model>::make(dict);
    auto badTagger        = TaggerFactory<hybrid_model>::make(dict);
    EXPECT_TRUE(badTagger == nullptr);
}


using Param   = PHARE::core::Span<double const>;
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

    using phare_types = PHARE::core::PHARE_Types<opts>;
    using Electromag  = phare_types::Hybrid::Electromag_t;
    using Ions        = phare_types::Hybrid::Ions_t;
    using Electrons   = phare_types::Hybrid::Electrons_t;
    using GridLayoutT
        = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;

    struct SinglePatchHybridModel
    {
        using gridlayout_type           = GridLayoutT;
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
        , tags(20 + layout.options.field_ghost_width)
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



template<std::size_t dim, std::size_t interp>
using HybridLayout = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;

constexpr PHARE::SimOpts mhdOpts(std::size_t dim)
{
    return PHARE::SimOpts{dim,
                          0,
                          0,
                          PHARE::MHDOpts::ReconstructionType::Linear,
                          PHARE::MHDOpts::SlopeLimiterType::VanLeer,
                          PHARE::MHDOpts::RiemannSolverType::Rusanov};
}

template<std::size_t dim>
using MHDLayout = PHARE::core::PHARE_Types<mhdOpts(dim)>::MHD::GridLayout_t;


template<std::size_t dim, std::size_t interp = 1>
struct MockHybridModel
{
    static constexpr std::string_view model_type_name = "HybridModel";
    using gridlayout_type                             = HybridLayout<dim, interp>;
    static constexpr auto dimension                   = dim;

    struct Electromag
    {
        PHARE::core::UsableVecField<dim> E;
        PHARE::core::UsableVecField<dim> B;
    };

    struct State
    {
        Electromag electromag;
    };

    explicit MockHybridModel(std::uint32_t cells = 20)
        : layout{TestGridLayout<gridlayout_type>::make(cells)}
        , state{{{"EM_E", layout, PHARE::core::HybridQuantity::Vector::E},
                 {"EM_B", layout, PHARE::core::HybridQuantity::Vector::B}}}
    {
    }

    gridlayout_type layout;
    State state;
};


template<std::size_t dim>
struct MockMHDModel
{
    static constexpr std::string_view model_type_name = "MHDModel";
    using gridlayout_type                             = MHDLayout<dim>;
    static constexpr auto dimension                   = dim;

    explicit MockMHDModel(std::uint32_t cells = 20)
        : layout{TestGridLayout<gridlayout_type>::make(cells)}
        , usable{layout, std::string{"MHD"}}
        , state{std::string{"MHD"}}
    {
        usable.set_on(state);
    }

    gridlayout_type layout;
    PHARE::core::UsableMHDState<dim> usable;
    PHARE::core::MHDState<PHARE::core::VecFieldMHD<dim>> state;
};


template<typename Model>
std::string resolveError(Model& model, std::string const& name)
{
    try
    {
        tagging::resolveQuantity(model, name);
    }
    catch (std::runtime_error const& e)
    {
        return e.what();
    }
    return {};
}

bool contains(std::string const& s, std::string const& sub)
{
    return s.find(sub) != std::string::npos;
}


TEST(TaggingQuantities, parseQuantitiesReadsList)
{
    PHARE::initializer::PHAREDict dict;
    dict["nbr_quantities"]  = 2;
    dict["Q0"]["name"]      = std::string{"B"};
    dict["Q0"]["threshold"] = 0.1;
    dict["Q1"]["name"]      = std::string{"rho"};
    dict["Q1"]["threshold"] = 0.4;
    auto const expected     = std::vector<tagging::QuantityThreshold>{{"B", 0.1}, {"rho", 0.4}};
    EXPECT_TRUE(tagging::parseQuantities(dict) == expected);
}

TEST(TaggingQuantities, parseQuantitiesFallsBackToBAtThreshold)
{
    PHARE::initializer::PHAREDict dict;
    dict["threshold"] = 0.2;
    EXPECT_TRUE(tagging::parseQuantities(dict)
                == (std::vector<tagging::QuantityThreshold>{{"B", 0.2}}));

    PHARE::initializer::PHAREDict empty;
    EXPECT_TRUE(tagging::parseQuantities(empty)
                == (std::vector<tagging::QuantityThreshold>{{"B", 0.1}}));
}

TEST(TaggingQuantities, hybridB)
{
    using PHARE::core::Component;
    MockHybridModel<2> model;
    auto& B = model.state.electromag.B;

    auto const comps = tagging::resolveQuantity(model, "B");
    ASSERT_EQ(comps.size(), 3u);
    EXPECT_EQ(comps[0].field(), &B.getComponent(Component::X));
    EXPECT_EQ(comps[1].field(), &B.getComponent(Component::Y));
    EXPECT_EQ(comps[2].field(), &B.getComponent(Component::Z));

    auto const by = tagging::resolveQuantity(model, "By");
    ASSERT_EQ(by.size(), 1u);
    EXPECT_EQ(by[0].field(), &B.getComponent(Component::Y));
}

TEST(TaggingQuantities, componentCentering)
{
    using Model = MockHybridModel<2>;
    Model model;
    auto const bx = tagging::resolveQuantity(model, "Bx");
    ASSERT_EQ(bx.size(), 1u);
    EXPECT_EQ(bx[0].centering(),
              Model::gridlayout_type::centering(PHARE::core::HybridQuantity::Scalar::Bx));
}

TEST(TaggingQuantities, mhdAllowed)
{
    using PHARE::core::Component;
    MockMHDModel<2> model;
    auto& s = model.state;

    for (std::string const name : {"rho", "density"})
    {
        auto const c = tagging::resolveQuantity(model, name);
        ASSERT_EQ(c.size(), 1u) << name;
        EXPECT_EQ(c[0].field(), &s.rho) << name;
    }
    for (std::string const name : {"rhoV", "momentum"})
    {
        auto const c = tagging::resolveQuantity(model, name);
        ASSERT_EQ(c.size(), 3u) << name;
        EXPECT_EQ(c[0].field(), &s.rhoV.getComponent(Component::X)) << name;
        EXPECT_EQ(c[1].field(), &s.rhoV.getComponent(Component::Y)) << name;
        EXPECT_EQ(c[2].field(), &s.rhoV.getComponent(Component::Z)) << name;
    }

    auto const rhoVy = tagging::resolveQuantity(model, "rhoVy");
    ASSERT_EQ(rhoVy.size(), 1u);
    EXPECT_EQ(rhoVy[0].field(), &s.rhoV.getComponent(Component::Y));

    auto const etot = tagging::resolveQuantity(model, "Etot");
    ASSERT_EQ(etot.size(), 1u);
    EXPECT_EQ(etot[0].field(), &s.Etot);

    auto const b = tagging::resolveQuantity(model, "B");
    ASSERT_EQ(b.size(), 3u);
    EXPECT_EQ(b[0].field(), &s.B.getComponent(Component::X));
    EXPECT_EQ(b[1].field(), &s.B.getComponent(Component::Y));
    EXPECT_EQ(b[2].field(), &s.B.getComponent(Component::Z));
}

TEST(TaggingQuantities, notAvailableForTaggingYet)
{
    MockHybridModel<2> hybrid;
    for (std::string const name :
         {"J", "Jx", "Jy", "Jz", "E", "Ex", "Ey", "Ez", "density", "massDensity", "chargeDensity",
          "velocity", "bulkVel", "momentum", "pressure"})
        EXPECT_TRUE(contains(resolveError(hybrid, name), "not available for tagging yet")) << name;

    MockMHDModel<2> mhd;
    for (std::string const name : {"J", "Jx", "Jy", "Jz", "E", "Ex", "Ey", "Ez", "V", "Vx", "Vy",
                                   "Vz", "velocity", "P", "pressure"})
        EXPECT_TRUE(contains(resolveError(mhd, name), "not available for tagging yet")) << name;
}

TEST(TaggingQuantities, unknownListsAllowed)
{
    MockHybridModel<2> hybrid;
    auto const hybridError = resolveError(hybrid, "rho");
    EXPECT_TRUE(contains(hybridError, "unknown")) << hybridError;
    EXPECT_TRUE(contains(hybridError, "Bz")) << hybridError;

    MockMHDModel<2> mhd;
    auto const mhdError = resolveError(mhd, "nonsense");
    EXPECT_TRUE(contains(mhdError, "unknown")) << mhdError;
    EXPECT_TRUE(contains(mhdError, "Etot")) << mhdError;
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int testResult = RUN_ALL_TESTS();
    return testResult;
}
