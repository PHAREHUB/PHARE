

#include "phare_solver.hpp"
#include "simulator/simulator_def.hpp"

#include "amr/tagging/tagger_factory.hpp"
#include "amr/tagging/tagging_quantities.hpp"
#include "amr/tagging/default_tagger_strategy.hpp"
#include "amr/tagging/component_wise_tagger_strategy.hpp"
#include "amr/tagging/lohner_tagger_strategy.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <map>
#include <span>
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
auto resolved(Model& model, std::string const& name)
{
    return tagging::resolveQuantity<Model>(name)(model);
}

template<typename Model>
std::string resolveError(std::string const& name)
{
    try
    {
        tagging::resolveQuantity<Model>(name);
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
    auto const expected     = std::map<std::string, double>{{"B", 0.1}, {"rho", 0.4}};
    EXPECT_TRUE(tagging::parseQuantities(dict) == expected);
}

TEST(TaggingQuantities, parseQuantitiesFallsBackToBAtThreshold)
{
    PHARE::initializer::PHAREDict dict;
    dict["threshold"] = 0.2;
    EXPECT_TRUE(tagging::parseQuantities(dict) == (std::map<std::string, double>{{"B", 0.2}}));

    PHARE::initializer::PHAREDict empty;
    EXPECT_TRUE(tagging::parseQuantities(empty) == (std::map<std::string, double>{{"B", 0.1}}));
}

TEST(TaggingQuantities, parseQuantitiesRejectsDuplicates)
{
    PHARE::initializer::PHAREDict dict;
    dict["nbr_quantities"]  = 2;
    dict["Q0"]["name"]      = std::string{"B"};
    dict["Q0"]["threshold"] = 0.1;
    dict["Q1"]["name"]      = std::string{"B"};
    dict["Q1"]["threshold"] = 0.4;
    EXPECT_THROW(tagging::parseQuantities(dict), std::runtime_error);
}

TEST(TaggingQuantities, hybridB)
{
    using PHARE::core::Component;
    MockHybridModel<2> model;
    auto& B = model.state.electromag.B;

    auto const comps = resolved(model, "B");
    ASSERT_EQ(comps.size(), 3u);
    EXPECT_EQ(comps[0].data(), B.getComponent(Component::X).data());
    EXPECT_EQ(comps[1].data(), B.getComponent(Component::Y).data());
    EXPECT_EQ(comps[2].data(), B.getComponent(Component::Z).data());

    auto const by = resolved(model, "By");
    ASSERT_EQ(by.size(), 1u);
    EXPECT_EQ(by[0].data(), B.getComponent(Component::Y).data());
}

TEST(TaggingQuantities, mhdAllowed)
{
    using PHARE::core::Component;
    MockMHDModel<2> model;
    auto& s = model.state;

    for (std::string const name : {"rho", "density"})
    {
        auto const c = resolved(model, name);
        ASSERT_EQ(c.size(), 1u) << name;
        EXPECT_EQ(c[0].data(), s.rho.data()) << name;
    }
    for (std::string const name : {"rhoV", "momentum"})
    {
        auto const c = resolved(model, name);
        ASSERT_EQ(c.size(), 3u) << name;
        EXPECT_EQ(c[0].data(), s.rhoV.getComponent(Component::X).data()) << name;
        EXPECT_EQ(c[1].data(), s.rhoV.getComponent(Component::Y).data()) << name;
        EXPECT_EQ(c[2].data(), s.rhoV.getComponent(Component::Z).data()) << name;
    }

    auto const rhoVy = resolved(model, "rhoVy");
    ASSERT_EQ(rhoVy.size(), 1u);
    EXPECT_EQ(rhoVy[0].data(), s.rhoV.getComponent(Component::Y).data());

    auto const etot = resolved(model, "Etot");
    ASSERT_EQ(etot.size(), 1u);
    EXPECT_EQ(etot[0].data(), s.Etot.data());

    auto const b = resolved(model, "B");
    ASSERT_EQ(b.size(), 3u);
    EXPECT_EQ(b[0].data(), s.B.getComponent(Component::X).data());
    EXPECT_EQ(b[1].data(), s.B.getComponent(Component::Y).data());
    EXPECT_EQ(b[2].data(), s.B.getComponent(Component::Z).data());
}

TEST(TaggingQuantities, notAvailableForTaggingYet)
{
    for (std::string const name :
         {"J", "Jx", "Jy", "Jz", "E", "Ex", "Ey", "Ez", "density", "massDensity", "chargeDensity",
          "velocity", "bulkVel", "momentum", "pressure"})
        EXPECT_TRUE(
            contains(resolveError<MockHybridModel<2>>(name), "not available for tagging yet"))
            << name;

    for (std::string const name : {"J", "Jx", "Jy", "Jz", "E", "Ex", "Ey", "Ez", "V", "Vx", "Vy",
                                   "Vz", "velocity", "P", "pressure"})
        EXPECT_TRUE(contains(resolveError<MockMHDModel<2>>(name), "not available for tagging yet"))
            << name;
}

TEST(TaggingQuantities, unknownListsAllowed)
{
    auto const hybridError = resolveError<MockHybridModel<2>>("rho");
    EXPECT_TRUE(contains(hybridError, "unknown")) << hybridError;
    EXPECT_TRUE(contains(hybridError, "Bz")) << hybridError;

    auto const mhdError = resolveError<MockMHDModel<2>>("nonsense");
    EXPECT_TRUE(contains(mhdError, "unknown")) << mhdError;
    EXPECT_TRUE(contains(mhdError, "Etot")) << mhdError;
}



template<typename Field, typename Fn>
void fillField(Field& field, Fn&& fn)
{
    auto constexpr dim = Field::dimension;
    auto const shape   = field.shape();
    if constexpr (dim == 1)
        for (std::uint32_t i = 0; i < shape[0]; ++i)
            field(i) = fn(std::array<std::uint32_t, 1>{i});
    if constexpr (dim == 2)
        for (std::uint32_t i = 0; i < shape[0]; ++i)
            for (std::uint32_t j = 0; j < shape[1]; ++j)
                field(i, j) = fn(std::array<std::uint32_t, 2>{i, j});
    if constexpr (dim == 3)
        for (std::uint32_t i = 0; i < shape[0]; ++i)
            for (std::uint32_t j = 0; j < shape[1]; ++j)
                for (std::uint32_t k = 0; k < shape[2]; ++k)
                    field(i, j, k) = fn(std::array<std::uint32_t, 3>{i, j, k});
}

template<typename Model>
std::vector<int> sentinelTags(Model const& model)
{
    return std::vector<int>(PHARE::core::product(model.layout.nbrCells()), -1);
}

template<typename Model>
std::vector<int> strategyTags(Model& model, PHARE::initializer::PHAREDict const& dict)
{
    auto tags = sentinelTags(model);
    DefaultTaggerStrategy<Model>{dict}.tag(model, model.layout, tags.data());
    return tags;
}

PHARE::initializer::PHAREDict quantitiesDict(std::map<std::string, double> const& qs)
{
    PHARE::initializer::PHAREDict dict;
    dict["nbr_quantities"] = static_cast<int>(qs.size());
    std::size_t i          = 0;
    for (auto const& [name, threshold] : qs)
    {
        auto const path         = "Q" + std::to_string(i++);
        dict[path]["name"]      = name;
        dict[path]["threshold"] = threshold;
    }
    return dict;
}

template<typename Model>
void fillMHDStep(Model& model, double low, double high)
{
    using PHARE::core::Component;
    auto const ghosts = Model::gridlayout_type::options.field_ghost_width;
    for (std::size_t c = 0; c < 3; ++c)
        fillField(model.state.B.getComponent(static_cast<Component>(c)),
                  [](auto const&) { return 1.0; });
    fillField(model.state.rho, [&](auto const& idx) { return idx[0] < ghosts + 10 ? low : high; });
}

template<typename Model>
bool tagged(Model const& model, std::vector<int> const& tags, std::uint32_t ix, std::uint32_t iy)
{
    auto const nx = model.layout.nbrCells()[0];
    return tags[ix + nx * iy] == 1;
}

template<typename Model>
void expectColumns(Model const& model, std::vector<int> const& tags,
                   std::vector<std::uint32_t> const& taggedColumns,
                   std::vector<std::uint32_t> const& untaggedColumns)
{
    auto const ny = model.layout.nbrCells()[1];
    for (std::uint32_t iy = 0; iy < ny; ++iy)
    {
        for (auto ix : taggedColumns)
            EXPECT_TRUE(tagged(model, tags, ix, iy)) << ix << " " << iy;
        for (auto ix : untaggedColumns)
            EXPECT_FALSE(tagged(model, tags, ix, iy)) << ix << " " << iy;
    }
}

TEST(DefaultStrategy, TagsScalarStep)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 1., 2.);
    auto const tags = strategyTags(model, quantitiesDict({{"rho", 0.1}}));
    expectColumns(model, tags, {8, 9}, {0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19});
}

TEST(DefaultStrategy, UnionOverQuantities)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 1., 2.);
    auto const onlyB = strategyTags(model, quantitiesDict({{"B", 0.1}}));
    EXPECT_EQ(std::count(onlyB.begin(), onlyB.end(), 1), 0);
    auto const both = strategyTags(model, quantitiesDict({{"B", 0.1}, {"rho", 0.1}}));
    EXPECT_NE(std::count(both.begin(), both.end(), 1), 0);
    EXPECT_EQ(both, strategyTags(model, quantitiesDict({{"rho", 0.1}})));
}

TEST(DefaultStrategy, PerQuantityThreshold)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 1., 2.);
    auto const low = strategyTags(model, quantitiesDict({{"rho", 0.1}}));
    EXPECT_NE(std::count(low.begin(), low.end(), 1), 0);
    auto const high = strategyTags(model, quantitiesDict({{"rho", 1e10}}));
    EXPECT_EQ(std::count(high.begin(), high.end(), 1), 0);
}



template<typename Model, std::uint32_t reach>
class ConstantCriterion : public ComponentWiseTaggerStrategy<ConstantCriterion<Model, reach>, Model>
{
    using Super      = ComponentWiseTaggerStrategy<ConstantCriterion<Model, reach>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = reach;

    ConstantCriterion(PHARE::initializer::PHAREDict const& dict, double value)
        : Super{dict}
        , value_{value}
    {
    }

    double indicator(field_type const&, std::array<std::uint32_t, Model::dimension> const&) const
    {
        return value_;
    }

private:
    double value_;
};

template<typename Model>
class CellValueCriterion : public ComponentWiseTaggerStrategy<CellValueCriterion<Model>, Model>
{
    using Super      = ComponentWiseTaggerStrategy<CellValueCriterion<Model>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = 0;

    explicit CellValueCriterion(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    double indicator(field_type const& field,
                     std::array<std::uint32_t, Model::dimension> const& cell) const
    {
        return field(cell);
    }
};

template<typename Model>
class NormCriterion : public ComponentWiseTaggerStrategy<NormCriterion<Model>, Model>
{
    using Super      = ComponentWiseTaggerStrategy<NormCriterion<Model>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = 0;

    explicit NormCriterion(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    double indicator(field_type const& field,
                     std::array<std::uint32_t, Model::dimension> const& cell) const
    {
        return field(cell);
    }

    static double combine(std::span<double const> indicators)
    {
        double sum = 0.;
        for (auto const c : indicators)
            sum += c * c;
        return std::sqrt(sum);
    }
};

static_assert(tagging::max_tagged_components_v<MockHybridModel<2>> == 3);
static_assert(tagging::max_tagged_components_v<MockMHDModel<2>> == 3);

template<typename Strategy, typename Model>
std::vector<int> tagsOf(Strategy const& strategy, Model& model)
{
    auto tags = sentinelTags(model);
    strategy.tag(model, model.layout, tags.data());
    return tags;
}

bool allEqual(std::vector<int> const& tags, int value)
{
    return std::all_of(tags.begin(), tags.end(), [&](int tag) { return tag == value; });
}

template<typename Field>
void fillConstant(Field& field, double value)
{
    fillField(field, [&](auto const&) { return value; });
}

template<typename Model>
void fillB(Model& model, double bx, double by, double bz)
{
    using PHARE::core::Component;
    auto& B = model.state.electromag.B;
    fillConstant(B.getComponent(Component::X), bx);
    fillConstant(B.getComponent(Component::Y), by);
    fillConstant(B.getComponent(Component::Z), bz);
}

template<typename Model>
auto constexpr ghostWidth
    = static_cast<std::uint32_t>(Model::gridlayout_type::options.field_ghost_width);


TEST(ComponentWise, OrOverComponentsOfAVector)
{
    using Model = MockHybridModel<2>;
    Model model;
    fillB(model, 0., 0., 1.);
    EXPECT_TRUE(
        allEqual(tagsOf(CellValueCriterion<Model>{quantitiesDict({{"B", 0.5}})}, model), 1));
    EXPECT_TRUE(
        allEqual(tagsOf(CellValueCriterion<Model>{quantitiesDict({{"Bx", 0.5}})}, model), 0));
}

TEST(ComponentWise, PerQuantityThresholds)
{
    using Model = MockMHDModel<2>;
    Model model;
    fillConstant(model.state.rho, 1.);
    fillConstant(model.state.Etot, 3.);
    auto const tags = [&](std::map<std::string, double> const& quantities) {
        return tagsOf(CellValueCriterion<Model>{quantitiesDict(quantities)}, model);
    };
    EXPECT_TRUE(allEqual(tags({{"rho", 2.}, {"Etot", 2.}}), 1));
    EXPECT_TRUE(allEqual(tags({{"rho", 2.}, {"Etot", 5.}}), 0));
    EXPECT_TRUE(allEqual(tags({{"rho", 1.}}), 0));
}

template<typename Model>
void fillMHDB(Model& model, double bx, double by, double bz)
{
    using PHARE::core::Component;
    fillConstant(model.state.B.getComponent(Component::X), bx);
    fillConstant(model.state.B.getComponent(Component::Y), by);
    fillConstant(model.state.B.getComponent(Component::Z), bz);
}

TEST(ComponentWise, CombineReducesTheComponentsOfOneQuantity)
{
    using Model = MockMHDModel<2>;
    Model model;
    auto const normTags = [&](std::map<std::string, double> const& quantities) {
        return tagsOf(NormCriterion<Model>{quantitiesDict(quantities)}, model);
    };

    fillMHDB(model, 0.4, 0.4, 0.);
    fillConstant(model.state.rho, 0.);
    EXPECT_TRUE(
        allEqual(tagsOf(CellValueCriterion<Model>{quantitiesDict({{"B", 0.5}})}, model), 0));
    EXPECT_TRUE(allEqual(normTags({{"B", 0.5}}), 1));

    fillMHDB(model, 0.3, 0.3, 0.);
    fillConstant(model.state.rho, 0.3);
    EXPECT_TRUE(allEqual(normTags({{"B", 0.5}, {"rho", 0.5}}), 0));
    EXPECT_TRUE(allEqual(normTags({{"rho", 0.2}}), 1));
}

TEST(ComponentWise, CombineNaNOnlyFailsItsQuantity)
{
    using Model = MockMHDModel<2>;
    Model model;
    fillMHDB(model, std::numeric_limits<double>::quiet_NaN(), 1., 0.);
    fillConstant(model.state.rho, 1.);
    EXPECT_TRUE(allEqual(tagsOf(NormCriterion<Model>{quantitiesDict({{"B", 0.5}})}, model), 0));
    EXPECT_TRUE(allEqual(
        tagsOf(NormCriterion<Model>{quantitiesDict({{"B", 0.5}, {"rho", 0.5}})}, model), 1));
}

TEST(ComponentWise, NaNComponentDoesNotMaskOthers)
{
    using Model = MockHybridModel<2>;
    Model model;
    fillB(model, std::numeric_limits<double>::quiet_NaN(), 0., 1.);
    EXPECT_TRUE(
        allEqual(tagsOf(CellValueCriterion<Model>{quantitiesDict({{"B", 0.5}})}, model), 1));
}

TEST(ComponentWise, ReusedAcrossPatches)
{
    using Model = MockHybridModel<2>;
    Model patchA, patchB;
    fillB(patchA, 0., 0., 1.);
    fillB(patchB, 0., 0., 0.);
    CellValueCriterion<Model> const criterion{quantitiesDict({{"B", 0.5}})};
    EXPECT_TRUE(allEqual(tagsOf(criterion, patchA), 1));
    EXPECT_TRUE(allEqual(tagsOf(criterion, patchB), 0));
}

TEST(ComponentWise, InvalidQuantityThrowsAtConstruction)
{
    using Model = MockMHDModel<2>;
    EXPECT_THROW(CellValueCriterion<Model>{quantitiesDict({{"rho", 0.1}, {"nonsense", 0.1}})},
                 std::runtime_error);
    EXPECT_THROW(DefaultTaggerStrategy<Model>{quantitiesDict({{"V", 0.1}})}, std::runtime_error);
}

TEST(ComponentWise, TagsFollowFortranOrdering)
{
    {
        using Model = MockMHDModel<2>;
        Model model;
        fillField(model.state.rho,
                  [](auto const& idx) { return idx[0] >= ghostWidth<Model> + 10 ? 1. : 0.; });
        auto const tags = tagsOf(CellValueCriterion<Model>{quantitiesDict({{"rho", 0.5}})}, model);
        expectColumns(model, tags, {10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    }
    {
        using Model = MockMHDModel<3>;
        Model model{12};
        fillField(model.state.rho,
                  [](auto const& idx) { return idx[0] >= ghostWidth<Model> + 6 ? 1. : 0.; });
        auto const tags = tagsOf(CellValueCriterion<Model>{quantitiesDict({{"rho", 0.5}})}, model);
        auto const n    = model.layout.nbrCells();
        for (std::uint32_t ix = 0; ix < n[0]; ++ix)
            for (std::uint32_t iy = 0; iy < n[1]; ++iy)
                for (std::uint32_t iz = 0; iz < n[2]; ++iz)
                    EXPECT_EQ(tags[ix + n[0] * (iy + n[1] * iz)], ix >= 6 ? 1 : 0)
                        << ix << " " << iy << " " << iz;
    }
}


template<typename Model>
struct ComponentWiseWideStencil : public ::testing::Test
{
};

using HybridModels = testing::Types<MockHybridModel<1>, MockHybridModel<2>, MockHybridModel<3>>;
TYPED_TEST_SUITE(ComponentWiseWideStencil, HybridModels);

TYPED_TEST(ComponentWiseWideStencil, LeavesBorderUntagged)
{
    using Model        = TypeParam;
    auto constexpr dim = Model::dimension;
    auto constexpr gw  = ghostWidth<Model>;
    Model model{8};
    auto const dict     = quantitiesDict({{"B", 0.5}});
    auto const nbrCells = model.layout.nbrCells();

    auto const tags = tagsOf(ConstantCriterion<Model, gw + 1>{dict, 1.}, model);
    for (std::size_t i = 0; i < tags.size(); ++i)
    {
        auto rest     = i;
        bool interior = true;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const t = rest % nbrCells[d];
            rest /= nbrCells[d];
            interior = interior and t >= 1 and t + 1 < nbrCells[d];
        }
        EXPECT_EQ(tags[i], interior ? 1 : 0) << i;
    }

    EXPECT_TRUE(allEqual(tagsOf(ConstantCriterion<Model, gw>{dict, 1.}, model), 1));
}

TEST(ComponentWise, PatchNarrowerThanBorderIsAllUntagged)
{
    using Model = MockHybridModel<1>;
    Model model{2};
    auto const tags = tagsOf(
        ConstantCriterion<Model, ghostWidth<Model> + 2>{quantitiesDict({{"B", 0.5}}), 1.}, model);
    EXPECT_TRUE(allEqual(tags, 0));
}

TEST(ComponentWise, WarnsAtConstructionWhenStencilTooWide)
{
    using Model       = MockHybridModel<2>;
    auto constexpr gw = ghostWidth<Model>;
    Model model;
    auto const dict = quantitiesDict({{"B", 0.5}});

    testing::internal::CaptureStdout();
    ConstantCriterion<Model, gw + 1> const wide{dict, 1.};
    EXPECT_TRUE(contains(testing::internal::GetCapturedStdout(), "stencil reach"));

    testing::internal::CaptureStdout();
    ConstantCriterion<Model, gw> const fits{dict, 1.};
    EXPECT_EQ(testing::internal::GetCapturedStdout(), "");

    testing::internal::CaptureStdout();
    auto const tags = tagsOf(wide, model);
    EXPECT_EQ(testing::internal::GetCapturedStdout(), "");
}



TEST(DefaultStrategy, OneDimensionalNormOverComponents)
{
    using PHARE::core::Component;
    using Model = MockHybridModel<1>;
    Model model;
    auto const step = [](auto const& idx) { return idx[0] >= ghostWidth<Model> + 10 ? 2. : 0.; };
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(Component::Y), step);
    fillField(model.state.electromag.B.getComponent(Component::Z), step);

    auto const tags = strategyTags(model, quantitiesDict({{"B", 0.5}}));
    for (std::size_t i = 0; i < tags.size(); ++i)
        EXPECT_EQ(tags[i], i == 7 ? 1 : 0) << i;
    EXPECT_TRUE(allEqual(strategyTags(model, quantitiesDict({{"By", 0.5}})), 0));
}

TEST(DefaultStrategy, ReadsCellCentredValues)
{
    using PHARE::core::Component;
    using Model = MockHybridModel<2>;
    Model model;
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(Component::X),
              [](auto const& idx) { return idx[0] == ghostWidth<Model> + 10 ? 0.3 : 0.; });
    EXPECT_TRUE(allEqual(strategyTags(model, quantitiesDict({{"Bx", 0.2}})), 0));
}

template<typename Model>
std::pair<std::vector<int>, std::string> alternatingByTagsAndWarning()
{
    Model model;
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
              [](auto const& idx) { return idx[0] % 2 ? 10. : -10.; });
    auto tags = sentinelTags(model);
    testing::internal::CaptureStdout();
    DefaultTaggerStrategy<Model> const strategy{quantitiesDict({{"By", 0.5}})};
    auto const warning = testing::internal::GetCapturedStdout();
    strategy.tag(model, model.layout, tags.data());
    return {tags, warning};
}

TEST(DefaultStrategy, OneDimensionalBorderDependsOnGhostWidth)
{
    auto const [twoGhosts, twoGhostsWarning] = alternatingByTagsAndWarning<MockHybridModel<1, 1>>();
    EXPECT_EQ(twoGhosts.front(), 0);
    EXPECT_EQ(twoGhosts.back(), 0);
    EXPECT_TRUE(allEqual(std::vector<int>(twoGhosts.begin() + 1, twoGhosts.end() - 1), 1));
    EXPECT_TRUE(contains(twoGhostsWarning, "stencil reach"));

    auto const [fourGhosts, fourGhostsWarning]
        = alternatingByTagsAndWarning<MockHybridModel<1, 2>>();
    EXPECT_TRUE(allEqual(fourGhosts, 1));
    EXPECT_EQ(fourGhostsWarning, "");
}

TEST(DefaultStrategy, TagsScalarStep3DAlongEachAxis)
{
    using Model = MockMHDModel<3>;
    Model model{12};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        fillField(model.state.rho,
                  [&](auto const& idx) { return idx[axis] >= ghostWidth<Model> + 6 ? 2. : 1.; });
        auto const tags = strategyTags(model, quantitiesDict({{"rho", 0.1}}));
        auto const n    = model.layout.nbrCells();
        for (std::uint32_t ix = 0; ix < n[0]; ++ix)
            for (std::uint32_t iy = 0; iy < n[1]; ++iy)
                for (std::uint32_t iz = 0; iz < n[2]; ++iz)
                {
                    auto const t = std::array{ix, iy, iz}[axis];
                    EXPECT_EQ(tags[ix + n[0] * (iy + n[1] * iz)], t == 4 or t == 5 ? 1 : 0)
                        << axis << ": " << ix << " " << iy << " " << iz;
                }
    }
}

TEST(DefaultStrategy, TwoDimensionalTestsComponentsSeparately)
{
    using PHARE::core::Component;
    using Model = MockHybridModel<2>;
    Model model;
    auto const step = [](auto const& idx) { return idx[0] >= ghostWidth<Model> + 10 ? 0.4 : 0.; };
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(Component::Y), step);
    fillField(model.state.electromag.B.getComponent(Component::Z), step);
    EXPECT_TRUE(allEqual(strategyTags(model, quantitiesDict({{"B", 0.5}})), 0));
    expectColumns(model, strategyTags(model, quantitiesDict({{"B", 0.3}})), {8},
                  {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
}


template<typename Fn>
double lohnerOnDualBy(Fn&& fn, std::uint32_t i, double reltol)
{
    using Model      = MockHybridModel<1>;
    using field_type = LohnerTaggerStrategy<Model>::recentered_field_type;
    Model model;
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    fillField(By, [&](auto const& idx) { return fn(idx[0]); });

    PHARE::initializer::PHAREDict dict;
    dict["params"]["reltol"] = reltol;
    dict["params"]["abstol"] = 1e-30;
    LohnerTaggerStrategy<Model> const lohner{dict};
    return lohner.indicator(field_type{By, {PHARE::core::QtyCentering::dual}}, {i});
}

TEST(LohnerIndicator, FlatAndLinearRampAreZero)
{
    EXPECT_EQ(lohnerOnDualBy([](auto) { return 3.0; }, 10, 0.02), 0.0);
    EXPECT_EQ(lohnerOnDualBy([](auto i) { return static_cast<double>(i); }, 10, 0.02), 0.0);
}

TEST(LohnerIndicator, StepNearOne)
{
    auto const step = [](auto i) { return i <= 10 ? 0.0 : 1.0; };
    EXPECT_DOUBLE_EQ(lohnerOnDualBy(step, 10, 0.02), 1.0 / 1.02);
}

TEST(LohnerIndicator, ReltolDampsSmallVariations)
{
    auto const step = [](auto i) { return i <= 10 ? 100.0 : 100.01; };
    EXPECT_LT(lohnerOnDualBy(step, 10, 0.02), 0.01);
    EXPECT_GT(lohnerOnDualBy(step, 10, 0.0), 0.9);
}


template<typename Model>
std::vector<int> lohnerTags(Model& model, PHARE::initializer::PHAREDict const& dict)
{
    auto tags = sentinelTags(model);
    LohnerTaggerStrategy<Model>{dict}.tag(model, model.layout, tags.data());
    return tags;
}

TEST(LohnerStrategy, TagsStepNotRamp)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 1., 2.);
    auto const step = lohnerTags(model, quantitiesDict({{"rho", 0.5}}));
    expectColumns(model, step, {9, 10}, {0, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 19});

    fillField(model.state.rho, [](auto const& idx) { return static_cast<double>(idx[0]); });
    auto const ramp = lohnerTags(model, quantitiesDict({{"rho", 0.5}}));
    EXPECT_EQ(std::count(ramp.begin(), ramp.end(), 1), 0);
    EXPECT_EQ(std::count(ramp.begin(), ramp.end(), -1), 0);
}

TEST(LohnerStrategy, ParamsFromDict)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 100., 100.01);
    auto const damped = lohnerTags(model, quantitiesDict({{"rho", 0.5}}));
    EXPECT_EQ(std::count(damped.begin(), damped.end(), 1), 0);

    auto dict                = quantitiesDict({{"rho", 0.5}});
    dict["params"]["reltol"] = 0.0;
    auto const undamped      = lohnerTags(model, dict);
    expectColumns(model, undamped, {9, 10}, {0, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 19});
}

TEST(LohnerStrategy, VectorQuantityIsOrOverComponents)
{
    using Model = MockHybridModel<1>;
    Model model;
    fillB(model, 1., 0., 0.);
    fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
              [](auto const& idx) { return idx[0] % 2 ? 1e-3 : -1e-3; });

    auto onB                 = quantitiesDict({{"B", 0.5}});
    onB["params"]["abstol"]  = 1e-8;
    auto onBx                = quantitiesDict({{"Bx", 0.5}});
    onBx["params"]["abstol"] = 1e-8;
    EXPECT_TRUE(allEqual(lohnerTags(model, onB), 1));
    EXPECT_TRUE(allEqual(lohnerTags(model, onBx), 0));
}

template<typename Model>
void fillAlternatingBy(Model& model, double amplitude)
{
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
              [&](auto const& idx) { return idx[0] % 2 ? amplitude : -amplitude; });
}

TEST(LohnerStrategy, AbstolIsANoiseFloor)
{
    MockHybridModel<1> model;
    fillAlternatingBy(model, 1e-3);

    auto below                = quantitiesDict({{"By", 0.5}});
    below["params"]["abstol"] = 1e-8;
    auto above                = quantitiesDict({{"By", 0.5}});
    above["params"]["abstol"] = 1e-2;
    EXPECT_TRUE(allEqual(lohnerTags(model, below), 1));
    EXPECT_TRUE(allEqual(lohnerTags(model, above), 0));
}

TEST(LohnerStrategy, AbstolMustBeSeveralTimesTheNoise)
{
    MockHybridModel<1> model;
    double constexpr sigma     = 1e-3;
    double constexpr threshold = 0.5;
    fillAlternatingBy(model, sigma);

    auto atNoise                     = quantitiesDict({{"By", threshold}});
    atNoise["params"]["abstol"]      = sigma;
    auto severalTimes                = quantitiesDict({{"By", threshold}});
    severalTimes["params"]["abstol"] = 4. * sigma * (1. / threshold - 1.);
    EXPECT_TRUE(allEqual(lohnerTags(model, atNoise), 1));
    EXPECT_TRUE(allEqual(lohnerTags(model, severalTimes), 0));
}

TEST(LohnerStrategy, DefaultAbstolIs1e8)
{
    MockHybridModel<1> model;
    fillAlternatingBy(model, 1e-9);

    auto noFloor                = quantitiesDict({{"By", 0.5}});
    noFloor["params"]["abstol"] = 0.;
    EXPECT_TRUE(allEqual(lohnerTags(model, quantitiesDict({{"By", 0.5}})), 0));
    EXPECT_TRUE(allEqual(lohnerTags(model, noFloor), 1));
}

TEST(TaggerFactory, hybridAndMHDMethods)
{
    using hybrid_model = PHARE::solver::PHARE_Types<PHARE::SimOpts{1ul, 1ul, 2ul}>::Hybrid::Model_t;
    using mhd_model    = PHARE::solver::PHARE_Types<mhdOpts(1)>::MHD::Model_t;

    for (std::string const method : {"default", "lohner"})
    {
        PHARE::initializer::PHAREDict dict;
        dict["hybrid_method"] = method;
        dict["mhd_method"]    = method;
        EXPECT_TRUE(TaggerFactory<hybrid_model>::make(dict) != nullptr) << method;
        EXPECT_TRUE(TaggerFactory<mhd_model>::make(dict) != nullptr) << method;
    }

    PHARE::initializer::PHAREDict dict;
    dict["hybrid_method"] = std::string{"nope"};
    dict["mhd_method"]    = std::string{"nope"};
    EXPECT_TRUE(TaggerFactory<hybrid_model>::make(dict) == nullptr);
    EXPECT_TRUE(TaggerFactory<mhd_model>::make(dict) == nullptr);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int testResult = RUN_ALL_TESTS();
    return testResult;
}
