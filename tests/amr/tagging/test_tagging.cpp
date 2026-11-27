

#include "initializer/data_provider.hpp"

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/models/mhd_state.hpp"
#include "core/models/options/mhd_options.hpp"
#include "core/models/options/mhd_options_def.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/utilities/array.hpp"
#include "core/utilities/numeric.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/span.hpp"
#include "phare_core.hpp"

#include "amr/tagging/default_tagger_strategy.hpp"
#include "amr/tagging/gradient_scale_tagger_strategy.hpp"
#include "amr/tagging/lohner_cross_derivatives_tagger_strategy.hpp"
#include "amr/tagging/lohner_parameters.hpp"
#include "amr/tagging/lohner_per_direction_tagger_strategy.hpp"
#include "amr/tagging/per_quantity_tagger_strategy.hpp"
#include "amr/tagging/tagger_factory.hpp"
#include "amr/tagging/tagging_quantities.hpp"
#include "amr/types/amr_types.hpp"
#include "phare_solver.hpp"

#include "simulator/simulator_def.hpp"

#include "dict.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

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


// minimal model for the tagger: only the magnetic field is read
struct UniformBModel
{
    using gridlayout_type
        = PHARE::core::PHARE_Types<PHARE::SimOpts{1ul, 1ul}>::Hybrid::GridLayout_t;
    static constexpr std::size_t dimension            = 1;
    static constexpr std::string_view model_type_name = "HybridModel";
    struct
    {
        struct
        {
            UsableVecField<dimension>& B;
        } electromag;
    } state;
};

TEST(test_tagger, everyCellIsWrittenIn1D)
{
    // SAMRAI pre-sets to 1 the tags of cells covered by the existing finer level before calling
    // the tagger, so every cell must be overwritten or it stays refined forever.
    // At interp order 1 the last cell of the patch is not evaluated: it must still be reset.
    auto static constexpr dim = UniformBModel::dimension;
    using GridLayoutT         = UniformBModel::gridlayout_type;

    std::uint32_t constexpr nCells = 20;
    auto layout                    = TestGridLayout<GridLayoutT>::make(nCells);
    UsableVecField<dim> B{"EM_B", layout, HybridQuantity::Vector::B}; // uniform (zero) field

    UniformBModel model{{{B}}};

    PHARE::initializer::PHAREDict dict;
    dict["threshold"] = 0.1;
    DefaultTaggerStrategy<UniformBModel> tagger{dict};

    std::vector<int> tags(nCells, 1); // as if a finer level covered the whole patch
    tagger.tag(model, layout, tags.data());

    for (std::size_t iCell = 0; iCell < nCells; ++iCell)
        EXPECT_EQ(tags[iCell], 0) << "cell " << iCell;
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
    dict["Q0"]["reltol"]    = 0.02;
    dict["Q1"]["name"]      = std::string{"rho"};
    dict["Q1"]["threshold"] = 0.4;
    auto const parsed       = tagging::parseQuantities(dict);
    ASSERT_EQ(parsed.size(), 2u);
    EXPECT_EQ(parsed.at("B").threshold, 0.1);
    EXPECT_EQ(parsed.at("rho").threshold, 0.4);
    EXPECT_EQ(cppdict::get_value(parsed.at("B").parameters, "reltol", 0.), 0.02);
    EXPECT_EQ(cppdict::get_value(parsed.at("rho").parameters, "reltol", 0.), 0.);
}

TEST(TaggingQuantities, parseQuantitiesFallsBackToBAtThreshold)
{
    PHARE::initializer::PHAREDict dict;
    dict["threshold"] = 0.2;
    auto const parsed = tagging::parseQuantities(dict);
    ASSERT_EQ(parsed.size(), 1u);
    EXPECT_EQ(parsed.at("B").threshold, 0.2);

    PHARE::initializer::PHAREDict empty;
    auto const fallback = tagging::parseQuantities(empty);
    ASSERT_EQ(fallback.size(), 1u);
    EXPECT_EQ(fallback.at("B").threshold, 0.1);
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

PHARE::initializer::PHAREDict
quantitiesDict(std::map<std::string, double> const& qs,
               std::map<std::string, std::map<std::string, double>> const& parameters = {})
{
    PHARE::initializer::PHAREDict dict;
    dict["nbr_quantities"] = static_cast<int>(qs.size());
    std::size_t i          = 0;
    for (auto const& [name, threshold] : qs)
    {
        auto const path         = "Q" + std::to_string(i++);
        dict[path]["name"]      = name;
        dict[path]["threshold"] = threshold;
        if (auto const found = parameters.find(name); found != parameters.end())
            for (auto const& [key, value] : found->second)
                dict[path][key] = value;
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
            EXPECT_EQ(tags[ix + model.layout.nbrCells()[0] * iy], 0) << ix << " " << iy;
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
class ConstantCriterion : public PerQuantityTaggerStrategy<ConstantCriterion<Model, reach>, Model>
{
    using Super      = PerQuantityTaggerStrategy<ConstantCriterion<Model, reach>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = reach;

    ConstantCriterion(PHARE::initializer::PHAREDict const& dict, double value)
        : Super{dict}
        , value_{value}
    {
    }

    double indicator(std::span<field_type const>,
                     std::array<std::uint32_t, Model::dimension> const&,
                     Model::gridlayout_type const&, NoQuantityParameters const&) const
    {
        return value_;
    }

private:
    double value_;
};

template<typename Model>
class CellValueCriterion : public PerQuantityTaggerStrategy<CellValueCriterion<Model>, Model>
{
    using Super      = PerQuantityTaggerStrategy<CellValueCriterion<Model>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = 0;

    explicit CellValueCriterion(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    double indicator(std::span<field_type const> components,
                     std::array<std::uint32_t, Model::dimension> const& cell,
                     Model::gridlayout_type const&, NoQuantityParameters const&) const
    {
        auto value = -std::numeric_limits<double>::infinity();
        for (auto const& component : components)
            value = std::max(value, component(cell));
        return value;
    }
};

template<typename Model>
class NormCriterion : public PerQuantityTaggerStrategy<NormCriterion<Model>, Model>
{
    using Super      = PerQuantityTaggerStrategy<NormCriterion<Model>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = 0;

    explicit NormCriterion(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    double indicator(std::span<field_type const> components,
                     std::array<std::uint32_t, Model::dimension> const& cell,
                     Model::gridlayout_type const&, NoQuantityParameters const&) const
    {
        double sum = 0.;
        for (auto const& component : components)
            sum += component(cell) * component(cell);
        return std::sqrt(sum);
    }
};

template<typename Model>
class ComponentsCriterion : public PerQuantityTaggerStrategy<ComponentsCriterion<Model>, Model>
{
    using Super      = PerQuantityTaggerStrategy<ComponentsCriterion<Model>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = 0;

    explicit ComponentsCriterion(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    double indicator(std::span<field_type const> components,
                     std::array<std::uint32_t, Model::dimension> const& cell,
                     Model::gridlayout_type const&, NoQuantityParameters const&) const
    {
        if (components.size() == 3)
            return components[0](cell) == 1. and components[1](cell) == 2.
                           and components[2](cell) == 3.
                       ? 1.
                       : 0.;
        return components.size() == 1 and components[0](cell) == 4. ? 1. : 0.;
    }
};

template<typename Model>
class MeshSizeCriterion : public PerQuantityTaggerStrategy<MeshSizeCriterion<Model>, Model>
{
    using Super      = PerQuantityTaggerStrategy<MeshSizeCriterion<Model>, Model>;
    using field_type = Super::recentered_field_type;

public:
    static constexpr std::uint32_t stencil_reach = 0;

    explicit MeshSizeCriterion(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    double indicator(std::span<field_type const>,
                     std::array<std::uint32_t, Model::dimension> const&,
                     Model::gridlayout_type const& layout, NoQuantityParameters const&) const
    {
        return layout.meshSize()[0];
    }
};

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


TEST(PerQuantity, PerQuantityThresholds)
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

TEST(PerQuantity, IndicatorReceivesEveryComponentInOrder)
{
    using Model = MockMHDModel<2>;
    Model model;
    fillMHDB(model, 1., 2., 3.);
    fillConstant(model.state.rho, 4.);
    auto const tags = [&](std::map<std::string, double> const& quantities) {
        return tagsOf(ComponentsCriterion<Model>{quantitiesDict(quantities)}, model);
    };
    EXPECT_TRUE(allEqual(tags({{"B", 0.5}}), 1));
    EXPECT_TRUE(allEqual(tags({{"rho", 0.5}}), 1));
    fillMHDB(model, 3., 2., 1.);
    EXPECT_TRUE(allEqual(tags({{"B", 0.5}}), 0));
}

TEST(PerQuantity, NaNOnlyFailsItsQuantity)
{
    using Model = MockMHDModel<2>;
    Model model;
    fillMHDB(model, std::numeric_limits<double>::quiet_NaN(), 1., 0.);
    fillConstant(model.state.rho, 1.);
    EXPECT_TRUE(allEqual(tagsOf(NormCriterion<Model>{quantitiesDict({{"B", 0.5}})}, model), 0));
    EXPECT_TRUE(allEqual(
        tagsOf(NormCriterion<Model>{quantitiesDict({{"B", 0.5}, {"rho", 0.5}})}, model), 1));
}

TEST(PerQuantity, ReusedAcrossPatches)
{
    using Model = MockHybridModel<2>;
    Model patchA, patchB;
    fillB(patchA, 0., 0., 1.);
    fillB(patchB, 0., 0., 0.);
    CellValueCriterion<Model> const criterion{quantitiesDict({{"B", 0.5}})};
    EXPECT_TRUE(allEqual(tagsOf(criterion, patchA), 1));
    EXPECT_TRUE(allEqual(tagsOf(criterion, patchB), 0));
}

TEST(PerQuantity, IndicatorReceivesThePatchLayout)
{
    using Model = MockMHDModel<2>;
    Model model;
    auto const dx = model.layout.meshSize()[0];
    EXPECT_TRUE(
        allEqual(tagsOf(MeshSizeCriterion<Model>{quantitiesDict({{"rho", 0.9 * dx}})}, model), 1));
    EXPECT_TRUE(
        allEqual(tagsOf(MeshSizeCriterion<Model>{quantitiesDict({{"rho", 1.1 * dx}})}, model), 0));
}

TEST(PerQuantity, InvalidQuantityThrowsAtConstruction)
{
    using Model = MockMHDModel<2>;
    EXPECT_THROW(CellValueCriterion<Model>{quantitiesDict({{"rho", 0.1}, {"nonsense", 0.1}})},
                 std::runtime_error);
    EXPECT_THROW(DefaultTaggerStrategy<Model>{quantitiesDict({{"V", 0.1}})}, std::runtime_error);
}

TEST(PerQuantity, TagsFollowFortranOrdering)
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
struct PerQuantityWideStencil : public ::testing::Test
{
};

using HybridModels = testing::Types<MockHybridModel<1>, MockHybridModel<2>, MockHybridModel<3>>;
TYPED_TEST_SUITE(PerQuantityWideStencil, HybridModels);

TYPED_TEST(PerQuantityWideStencil, LeavesBorderUntagged)
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

TEST(PerQuantity, PatchNarrowerThanBorderIsAllUntagged)
{
    using Model = MockHybridModel<1>;
    Model model{2};
    auto const tags = tagsOf(
        ConstantCriterion<Model, ghostWidth<Model> + 2>{quantitiesDict({{"B", 0.5}}), 1.}, model);
    EXPECT_TRUE(allEqual(tags, 0));
}

TEST(PerQuantity, WarnsAtConstructionWhenStencilTooWide)
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

TEST(DefaultStrategy, NaNComponentDoesNotMaskOthers)
{
    using PHARE::core::Component;
    using Model = MockHybridModel<2>;
    Model model;
    auto const step = [](auto const& idx) { return idx[0] >= ghostWidth<Model> + 10 ? 0.4 : 0.; };
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(Component::Y), step);
    auto const finite = strategyTags(model, quantitiesDict({{"B", 0.3}}));
    fillConstant(model.state.electromag.B.getComponent(Component::Z),
                 std::numeric_limits<double>::quiet_NaN());
    auto const withNaN = strategyTags(model, quantitiesDict({{"B", 0.3}}));
    EXPECT_NE(std::count(finite.begin(), finite.end(), 1), 0);
    EXPECT_EQ(withNaN, finite);
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
    using field_type = LohnerPerDirectionTaggerStrategy<Model>::recentered_field_type;
    Model model;
    auto& By = model.state.electromag.B.getComponent(PHARE::core::Component::Y);
    fillField(By, [&](auto const& idx) { return fn(idx[0]); });

    LohnerPerDirectionTaggerStrategy<Model> const lohner{PHARE::initializer::PHAREDict{}};
    return lohner.indicator(
        std::array{field_type{By, {PHARE::core::QtyCentering::dual}, model.layout}}, {i},
        model.layout, LohnerParameters{reltol, 1e-30});
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
    LohnerPerDirectionTaggerStrategy<Model>{dict}.tag(model, model.layout, tags.data());
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

TEST(LohnerStrategy, ReltolFromQuantityDict)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 100., 100.01);
    auto const damped = lohnerTags(model, quantitiesDict({{"rho", 0.5}}));
    EXPECT_EQ(std::count(damped.begin(), damped.end(), 1), 0);

    auto const undamped
        = lohnerTags(model, quantitiesDict({{"rho", 0.5}}, {{"rho", {{"reltol", 0.0}}}}));
    expectColumns(model, undamped, {9, 10}, {0, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 19});
}

TEST(LohnerParameters, ReadsKeysWithDefaults)
{
    PHARE::initializer::PHAREDict q;
    q["reltol"]        = 0.02;
    q["ref_magnitude"] = 1e-6;
    auto const read    = LohnerParameters::read(q, "B");
    EXPECT_EQ(read.reltol, 0.02);
    EXPECT_EQ(read.refMagnitude, 1e-6);
    auto const defaults = LohnerParameters::read(PHARE::initializer::PHAREDict{}, "B");
    EXPECT_EQ(defaults.reltol, 0.01);
    EXPECT_EQ(defaults.refMagnitude, 1e-8);
}

std::string lohnerParametersError(std::string const& key, double value)
{
    PHARE::initializer::PHAREDict q;
    q[key] = value;
    try
    {
        LohnerParameters::read(q, "rho");
    }
    catch (std::runtime_error const& e)
    {
        return e.what();
    }
    return {};
}

TEST(LohnerParameters, RejectsNegativeValuesNamingTheQuantity)
{
    auto const reltol = lohnerParametersError("reltol", -0.1);
    EXPECT_TRUE(contains(reltol, "rho") and contains(reltol, "reltol")) << reltol;
    auto const refMagnitude = lohnerParametersError("ref_magnitude", -1.);
    EXPECT_TRUE(contains(refMagnitude, "rho") and contains(refMagnitude, "ref_magnitude"))
        << refMagnitude;
    EXPECT_TRUE(lohnerParametersError("reltol", 0.).empty());
    EXPECT_TRUE(lohnerParametersError("ref_magnitude", 0.).empty());
}

TEST(LohnerParameters, RejectsNonFiniteValues)
{
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();
    for (auto const& key : {std::string{"reltol"}, std::string{"ref_magnitude"}})
        for (double const value : {nan, inf})
        {
            auto const error = lohnerParametersError(key, value);
            EXPECT_TRUE(contains(error, "rho") and contains(error, key)) << key << " " << value;
        }
}

TEST(LohnerStrategy, PerQuantityParametersFollowTheName)
{
    MockMHDModel<2> model;
    fillMHDStep(model, 100., 100.01);
    fillConstant(model.state.Etot, 3.);
    PHARE::initializer::PHAREDict dict;
    dict["nbr_quantities"]  = 2;
    dict["Q0"]["name"]      = std::string{"rho"};
    dict["Q0"]["threshold"] = 0.5;
    dict["Q0"]["reltol"]    = 0.0;
    dict["Q1"]["name"]      = std::string{"Etot"};
    dict["Q1"]["threshold"] = 0.5;
    expectColumns(model, lohnerTags(model, dict), {9, 10},
                  {0, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 19});
}

template<typename Model>
void fillAlternatingBy(Model& model, double amplitude)
{
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
              [&](auto const& idx) { return idx[0] % 2 ? amplitude : -amplitude; });
}

TEST(LohnerStrategy, RefMagnitudeIsANoiseLevel)
{
    MockHybridModel<1> model;
    fillAlternatingBy(model, 1e-3);

    auto const below = quantitiesDict({{"By", 0.5}}, {{"By", {{"ref_magnitude", 1e-8}}}});
    auto const above = quantitiesDict({{"By", 0.5}}, {{"By", {{"ref_magnitude", 1e-2}}}});
    EXPECT_TRUE(allEqual(lohnerTags(model, below), 1));
    EXPECT_TRUE(allEqual(lohnerTags(model, above), 0));
}

TEST(LohnerStrategy, RefMagnitudeMustBeSeveralTimesTheNoise)
{
    MockHybridModel<1> model;
    double constexpr sigma     = 1e-3;
    double constexpr threshold = 0.5;
    fillAlternatingBy(model, sigma);

    auto const atNoise = quantitiesDict({{"By", threshold}}, {{"By", {{"ref_magnitude", sigma}}}});
    auto const severalTimes = quantitiesDict(
        {{"By", threshold}}, {{"By", {{"ref_magnitude", 4. * sigma * (1. / threshold - 1.)}}}});
    EXPECT_TRUE(allEqual(lohnerTags(model, atNoise), 1));
    EXPECT_TRUE(allEqual(lohnerTags(model, severalTimes), 0));
}

TEST(LohnerStrategy, DefaultRefMagnitudeIs1e8)
{
    MockHybridModel<1> model;
    fillAlternatingBy(model, 1e-9);

    auto const noRefMagnitude = quantitiesDict({{"By", 0.5}}, {{"By", {{"ref_magnitude", 0.}}}});
    EXPECT_TRUE(allEqual(lohnerTags(model, quantitiesDict({{"By", 0.5}})), 0));
    EXPECT_TRUE(allEqual(lohnerTags(model, noRefMagnitude), 1));
}

template<typename Strategy, typename Model, typename Fn>
double indicatorAt(Model& model, Fn&& fn, std::array<std::uint32_t, Model::dimension> cell,
                   double reltol, double refMagnitude,
                   typename Model::gridlayout_type const& layout)
{
    fillField(model.state.rho, [&](auto const& idx) { return fn(idx); });
    auto const allDual = PHARE::core::ConstArray<PHARE::core::QtyCentering, Model::dimension>(
        PHARE::core::QtyCentering::dual);
    return Strategy{PHARE::initializer::PHAREDict{}}.indicator(
        std::array{typename Strategy::recentered_field_type{model.state.rho, allDual, layout}},
        cell, layout, LohnerParameters{reltol, refMagnitude});
}

TEST(LohnerPerDirection, RatioOfSummedSquaresIn2D)
{
    MockMHDModel<2> model;
    auto const u = [](auto const& idx) {
        return static_cast<double>(idx[0] * idx[0] + 2 * idx[1] * idx[1]);
    };
    EXPECT_NEAR((indicatorAt<LohnerPerDirectionTaggerStrategy<MockMHDModel<2>>>(
                    model, u, {10, 10}, 0.01, 0., model.layout)),
                std::sqrt(20. / (52.02 * 52.02 + 92.04 * 92.04)), 1e-14);
}

TEST(LohnerPerDirection, ZeroFieldWithoutRefMagnitudeIsZero)
{
    MockMHDModel<2> model;
    EXPECT_EQ((indicatorAt<LohnerPerDirectionTaggerStrategy<MockMHDModel<2>>>(
                  model, [](auto const&) { return 0.; }, {10, 10}, 0.01, 0., model.layout)),
              0.);
}

TEST(LohnerPerDirection, DefaultReltolIs001)
{
    MockMHDModel<1> model;
    fillField(model.state.rho, [](auto const& idx) { return idx[0] <= 10 ? 0. : 1.; });
    auto const allDual
        = PHARE::core::ConstArray<PHARE::core::QtyCentering, 1>(PHARE::core::QtyCentering::dual);
    LohnerPerDirectionTaggerStrategy<MockMHDModel<1>> const lohner{PHARE::initializer::PHAREDict{}};
    using field_type = LohnerPerDirectionTaggerStrategy<MockMHDModel<1>>::recentered_field_type;
    EXPECT_NEAR(lohner.indicator(std::array{field_type{model.state.rho, allDual, model.layout}},
                                 {10}, model.layout, LohnerParameters{}),
                1. / 1.01, 1e-7);
}

using PerDirection2 = LohnerPerDirectionTaggerStrategy<MockMHDModel<2>>;
using Cross1        = LohnerCrossDerivativesTaggerStrategy<MockMHDModel<1>>;
using Cross2        = LohnerCrossDerivativesTaggerStrategy<MockMHDModel<2>>;
using Cross3        = LohnerCrossDerivativesTaggerStrategy<MockMHDModel<3>>;

TEST(LohnerCrossDerivatives, OneDimensionIsFlashEq83)
{
    MockMHDModel<1> model;
    auto const u = [](auto const& idx) { return static_cast<double>(idx[0] * idx[0]); };
    EXPECT_NEAR((indicatorAt<Cross1>(model, u, {10}, 0.01, 0., model.layout)), 8. / 84.08, 1e-14);
}

TEST(LohnerCrossDerivatives, SaddleSeenOnlyByCrossTerms)
{
    MockMHDModel<2> model;
    auto const saddle = [](auto const& idx) {
        return (static_cast<double>(idx[0]) - 10.) * (static_cast<double>(idx[1]) - 10.);
    };
    EXPECT_EQ((indicatorAt<PerDirection2>(model, saddle, {10, 10}, 0.01, 1e-8, model.layout)), 0.);
    EXPECT_NEAR((indicatorAt<Cross2>(model, saddle, {10, 10}, 0.01, 0., model.layout)), 4. / 4.04,
                1e-14);
}

TEST(LohnerCrossDerivatives, SaddleIn3DMatches2D)
{
    MockMHDModel<3> model;
    auto const saddle = [](auto const& idx) {
        return (static_cast<double>(idx[0]) - 10.) * (static_cast<double>(idx[1]) - 10.);
    };
    EXPECT_NEAR((indicatorAt<Cross3>(model, saddle, {10, 10, 10}, 0.01, 0., model.layout)),
                4. / 4.04, 1e-14);
}

TEST(LohnerCrossDerivatives, MeshSizeWeightsAreLevelIndependent)
{
    using Layout = MockMHDModel<2>::gridlayout_type;
    MockMHDModel<2> model;
    auto const u = [](auto const& idx) {
        auto const i = static_cast<double>(idx[0]), j = static_cast<double>(idx[1]);
        return i * i + (i - 10.) * (j - 10.);
    };
    auto const at = [&](double dx, double dy) {
        Layout const layout{{dx, dy}, {20u, 20u}, PHARE::core::Point<double, 2>{0., 0.}};
        return indicatorAt<Cross2>(model, u, {10, 10}, 0.01, 0., layout);
    };
    EXPECT_NEAR(at(1., 1.), 0.08218558436423791, 1e-14);
    EXPECT_NEAR(at(1., 2.), 0.09018587026447226, 1e-14);
    EXPECT_NEAR(at(2., 4.), at(1., 2.), 1e-14);
}

TEST(LohnerCrossDerivatives, ZeroFieldWithoutRefMagnitudeIsZero)
{
    MockMHDModel<2> model;
    EXPECT_EQ((indicatorAt<Cross2>(
                  model, [](auto const&) { return 0.; }, {10, 10}, 0.01, 0., model.layout)),
              0.);
}

TEST(LohnerCrossDerivatives, DefaultReltolIs001)
{
    MockMHDModel<1> model;
    fillField(model.state.rho, [](auto const& idx) { return idx[0] <= 11 ? 0. : 1.; });
    auto const allDual
        = PHARE::core::ConstArray<PHARE::core::QtyCentering, 1>(PHARE::core::QtyCentering::dual);
    Cross1 const cross{PHARE::initializer::PHAREDict{}};
    EXPECT_NEAR(cross.indicator(std::array{Cross1::recentered_field_type{model.state.rho, allDual,
                                                                         model.layout}},
                                {10}, model.layout, LohnerParameters{}),
                1. / 1.01, 1e-7);
}

TEST(LohnerCrossDerivatives, BoundedOnRandomFields)
{
    MockMHDModel<2> model;
    std::mt19937 gen{42};
    std::uniform_real_distribution<double> dist{-5., 5.};
    fillField(model.state.rho, [&](auto const&) { return dist(gen); });
    auto const allDual
        = PHARE::core::ConstArray<PHARE::core::QtyCentering, 2>(PHARE::core::QtyCentering::dual);
    PHARE::initializer::PHAREDict dict;
    Cross2 const cross{dict};
    PerDirection2 const perDirection{dict};
    for (std::uint32_t i = 4; i < 20; ++i)
        for (std::uint32_t j = 4; j < 20; ++j)
        {
            EXPECT_LT(cross.indicator(std::array{Cross2::recentered_field_type{
                                          model.state.rho, allDual, model.layout}},
                                      {i, j}, model.layout, LohnerParameters{}),
                      1. / 1.01);
            EXPECT_LT(perDirection.indicator(std::array{PerDirection2::recentered_field_type{
                                                 model.state.rho, allDual, model.layout}},
                                             {i, j}, model.layout, LohnerParameters{}),
                      1. / 1.01);
        }
}

TEST(LohnerCrossDerivatives, CrossDerivativesTagsEveryCellAtGhostWidthTwo)
{
    using Model = MockHybridModel<2>;
    static_assert(ghostWidth<Model> == 2);
    Model model;
    auto const expectTagged = [&](std::uint32_t jump, std::vector<std::uint32_t> const& tagged) {
        fillB(model, 1., 0., 0.);
        fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
                  [&](auto const& idx) { return idx[0] >= ghostWidth<Model> + jump ? 1. : 0.; });
        auto tags = sentinelTags(model);
        LohnerCrossDerivativesTaggerStrategy<Model>{quantitiesDict({{"B", 0.5}})}.tag(
            model, model.layout, tags.data());
        std::vector<std::uint32_t> untagged;
        for (std::uint32_t ix = 0; ix < model.layout.nbrCells()[0]; ++ix)
            if (std::find(tagged.begin(), tagged.end(), ix) == tagged.end())
                untagged.push_back(ix);
        expectColumns(model, tags, tagged, untagged);
    };
    expectTagged(10, {8, 11});
    expectTagged(2, {0, 3});
    expectTagged(18, {16, 19});
}

template<typename Strategy, typename Model>
std::vector<int> methodTags(Model& model, PHARE::initializer::PHAREDict const& dict)
{
    auto tags = sentinelTags(model);
    Strategy{dict}.tag(model, model.layout, tags.data());
    return tags;
}

template<typename Model, typename Fx, typename Fy, typename Fz>
void fillRhoV(Model& model, Fx&& fx, Fy&& fy, Fz&& fz)
{
    using PHARE::core::Component;
    fillField(model.state.rhoV.getComponent(Component::X), fx);
    fillField(model.state.rhoV.getComponent(Component::Y), fy);
    fillField(model.state.rhoV.getComponent(Component::Z), fz);
}

template<typename Strategy, typename VecField, typename Layout, std::size_t dim,
         typename Parameters>
double vectorIndicatorAt(VecField const& vecfield, std::array<std::uint32_t, dim> const& cell,
                         Parameters const& parameters, Layout const& layout)
{
    using PHARE::core::Component;
    using field_type = Strategy::recentered_field_type;
    auto const allDual
        = PHARE::core::ConstArray<PHARE::core::QtyCentering, dim>(PHARE::core::QtyCentering::dual);
    std::array const components{field_type{vecfield.getComponent(Component::X), allDual, layout},
                                field_type{vecfield.getComponent(Component::Y), allDual, layout},
                                field_type{vecfield.getComponent(Component::Z), allDual, layout}};
    return Strategy{quantitiesDict({{"B", 1.}}, {{"B", {{"ref_magnitude", 1.}}}})}.indicator(
        components, cell, layout, parameters);
}

template<typename Model, typename Fn>
void forBothLohner(Fn&& fn)
{
    PHARE::initializer::PHAREDict const dict;
    fn(LohnerPerDirectionTaggerStrategy<Model>{dict});
    fn(LohnerCrossDerivativesTaggerStrategy<Model>{dict});
}

TEST(LohnerStrategy, SmallComponentOnStrongFieldIsNotTagged)
{
    using Model = MockHybridModel<1>;
    Model model;
    fillB(model, 1., 0., 0.);
    fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
              [](auto const& idx) { return (idx[0] / 2) % 2 ? -1e-3 : 1e-3; });
    forBothLohner<Model>([&](auto method) {
        using Strategy = decltype(method);
        EXPECT_TRUE(allEqual(methodTags<Strategy>(model, quantitiesDict({{"B", 0.5}})), 0));
        EXPECT_TRUE(allEqual(methodTags<Strategy>(model, quantitiesDict({{"By", 0.5}})), 1));
    });
}

TEST(LohnerStrategy, RotationAtConstantMagnitudeIsTagged)
{
    using PHARE::core::Component;
    using Model = MockHybridModel<1>;
    Model model;
    auto constexpr jump = ghostWidth<Model> + 10;
    fillB(model, 0., 0., 0.);
    fillField(model.state.electromag.B.getComponent(Component::Y),
              [](auto const& idx) { return idx[0] >= jump ? 0. : 1.; });
    fillField(model.state.electromag.B.getComponent(Component::Z),
              [](auto const& idx) { return idx[0] >= jump ? 1. : 0.; });
    auto const expected = std::sqrt(2.) / (std::sqrt(2.) + 4. * 0.01);

    auto const expectTagged = [&](auto method, std::vector<std::uint32_t> const& tagged) {
        using Strategy  = decltype(method);
        auto const tags = methodTags<Strategy>(model, quantitiesDict({{"B", 0.5}}));
        for (std::uint32_t i = 0; i < tags.size(); ++i)
        {
            auto const isTagged = std::find(tagged.begin(), tagged.end(), i) != tagged.end();
            EXPECT_EQ(tags[i], isTagged ? 1 : 0) << i;
        }
        for (auto const i : tagged)
            EXPECT_NEAR((vectorIndicatorAt<Strategy>(model.state.electromag.B,
                                                     std::array{i + ghostWidth<Model>},
                                                     LohnerParameters{}, model.layout)),
                        expected, 1e-7)
                << i;
    };
    expectTagged(LohnerPerDirectionTaggerStrategy<Model>{PHARE::initializer::PHAREDict{}}, {9, 10});
    expectTagged(LohnerCrossDerivativesTaggerStrategy<Model>{PHARE::initializer::PHAREDict{}},
                 {8, 9, 10, 11});
}

TEST(LohnerVector, PerDirectionHandComputed2D)
{
    MockMHDModel<2> model;
    fillRhoV(
        model, [](auto const& idx) { return static_cast<double>(idx[0] * idx[0]); },
        [](auto const& idx) { return static_cast<double>(2 * idx[1] * idx[1]); },
        [](auto const&) { return 0.; });
    auto const Dx = 40. + 0.01 * (std::sqrt(54641.) + 2. * std::sqrt(50000.) + std::sqrt(46561.));
    auto const Dy = 80. + 0.01 * (std::sqrt(68564.) + 2. * std::sqrt(50000.) + std::sqrt(36244.));
    EXPECT_NEAR((vectorIndicatorAt<PerDirection2>(model.state.rhoV, std::array{10u, 10u},
                                                  LohnerParameters{0.01, 0.}, model.layout)),
                std::sqrt(20. / (Dx * Dx + Dy * Dy)), 1e-14);
}

TEST(LohnerVector, CrossDerivativesHandComputed2D)
{
    MockMHDModel<2> model;
    fillRhoV(
        model, [](auto const& idx) { return static_cast<double>(idx[0] * idx[0]); },
        [](auto const& idx) { return static_cast<double>(2 * idx[1] * idx[1]); },
        [](auto const&) { return 0.; });
    EXPECT_NEAR((vectorIndicatorAt<Cross2>(model.state.rhoV, std::array{10u, 10u},
                                           LohnerParameters{0.01, 0.}, model.layout)),
                0.06617835648526077, 1e-14);
}

template<typename Strategy, typename Model, typename Fn>
void expectParallelMatchesScalar(Model& model, Fn&& f,
                                 std::array<std::uint32_t, Model::dimension> const& cell,
                                 typename Model::gridlayout_type const& layout)
{
    auto const scalar = indicatorAt<Strategy>(model, f, cell, 0.01, 0., layout);
    fillRhoV(
        model, [&](auto const& idx) { return 0.6 * f(idx); }, [](auto const&) { return 0.; },
        [&](auto const& idx) { return 0.8 * f(idx); });
    EXPECT_NEAR(
        (vectorIndicatorAt<Strategy>(model.state.rhoV, cell, LohnerParameters{0.01, 0.}, layout)),
        scalar, 1e-12);
}

TEST(LohnerVector, ParallelVectorMatchesScalar)
{
    {
        MockMHDModel<2> model;
        expectParallelMatchesScalar<PerDirection2>(
            model,
            [](auto const& idx) {
                return static_cast<double>(idx[0] * idx[0] + 2 * idx[1] * idx[1]);
            },
            {10, 10}, model.layout);
    }
    {
        using Layout = MockMHDModel<2>::gridlayout_type;
        MockMHDModel<2> model;
        Layout const layout{{1., 2.}, {20u, 20u}, PHARE::core::Point<double, 2>{0., 0.}};
        auto const u = [](auto const& idx) {
            auto const i = static_cast<double>(idx[0]), j = static_cast<double>(idx[1]);
            return i * i + (i - 10.) * (j - 10.);
        };
        expectParallelMatchesScalar<Cross2>(model, u, {10, 10}, layout);
        EXPECT_NEAR((vectorIndicatorAt<Cross2>(model.state.rhoV, std::array{10u, 10u},
                                               LohnerParameters{0.01, 0.}, layout)),
                    0.09018587026447226, 1e-12);
    }
    {
        MockMHDModel<3> model;
        auto const saddle = [](auto const& idx) {
            return (static_cast<double>(idx[0]) - 10.) * (static_cast<double>(idx[1]) - 10.);
        };
        expectParallelMatchesScalar<Cross3>(model, saddle, {10, 10, 10}, model.layout);
        EXPECT_NEAR((vectorIndicatorAt<Cross3>(model.state.rhoV, std::array{10u, 10u, 10u},
                                               LohnerParameters{0.01, 0.}, model.layout)),
                    4. / 4.04, 1e-12);
    }
}

TEST(LohnerVector, BoundedOnRandomVectorFields)
{
    MockMHDModel<2> model;
    std::mt19937 gen{42};
    std::uniform_real_distribution<double> dist{-5., 5.};
    auto const random = [&](auto const&) { return dist(gen); };
    fillRhoV(model, random, random, random);
    LohnerParameters const parameters;
    for (std::uint32_t i = 4; i < 20; ++i)
        for (std::uint32_t j = 4; j < 20; ++j)
        {
            EXPECT_LT((vectorIndicatorAt<PerDirection2>(model.state.rhoV, std::array{i, j},
                                                        parameters, model.layout)),
                      1. / 1.01);
            EXPECT_LT((vectorIndicatorAt<Cross2>(model.state.rhoV, std::array{i, j}, parameters,
                                                 model.layout)),
                      1. / 1.01);
        }
}

TEST(LohnerStrategy, NaNComponentFailsOnlyItsQuantity)
{
    using PHARE::core::Component;
    using Model = MockMHDModel<2>;
    Model model;
    fillMHDStep(model, 1., 2.);
    fillField(model.state.B.getComponent(Component::Y),
              [](auto const& idx) { return idx[0] < ghostWidth<Model> + 10 ? 1. : 2.; });
    fillConstant(model.state.B.getComponent(Component::X),
                 std::numeric_limits<double>::quiet_NaN());
    forBothLohner<Model>([&](auto method) {
        using Strategy   = decltype(method);
        auto const onB   = methodTags<Strategy>(model, quantitiesDict({{"B", 0.5}}));
        auto const onRho = methodTags<Strategy>(model, quantitiesDict({{"rho", 0.5}}));
        auto const both  = methodTags<Strategy>(model, quantitiesDict({{"B", 0.5}, {"rho", 0.5}}));
        EXPECT_TRUE(allEqual(onB, 0));
        EXPECT_NE(std::count(onRho.begin(), onRho.end(), 1), 0);
        EXPECT_EQ(both, onRho);
    });
}

TEST(LohnerStrategy, RoundoffVectorIsNotTagged)
{
    using Model = MockMHDModel<2>;
    Model model;
    std::mt19937 gen{7};
    std::uniform_real_distribution<double> dist{-1e-14, 1e-14};
    auto const random = [&](auto const&) { return dist(gen); };
    fillRhoV(model, random, random, random);
    forBothLohner<Model>([&](auto method) {
        using Strategy = decltype(method);
        EXPECT_TRUE(allEqual(methodTags<Strategy>(model, quantitiesDict({{"rhoV", 0.1}})), 0));
    });
}

template<typename Model, typename Fn>
double gradientScaleOnRho(Model& model, Fn&& fn, std::array<std::uint32_t, Model::dimension> cell,
                          double refMagnitude, typename Model::gridlayout_type const& layout)
{
    using Strategy = GradientScaleTaggerStrategy<Model>;
    fillField(model.state.rho, fn);
    auto const allDual = PHARE::core::ConstArray<PHARE::core::QtyCentering, Model::dimension>(
        PHARE::core::QtyCentering::dual);
    Strategy const strategy{
        quantitiesDict({{"rho", 1.}}, {{"rho", {{"ref_magnitude", refMagnitude}}}})};
    return strategy.indicator(
        std::array{typename Strategy::recentered_field_type{model.state.rho, allDual, layout}},
        cell, layout, GradientScaleParameters{refMagnitude});
}

TEST(GradientScale, ScalarRampIsHalfDifferenceOverValue)
{
    MockMHDModel<2> model;
    auto const ramp = [](auto const& idx) { return 3. + idx[0] + 2. * idx[1]; };
    EXPECT_NEAR(gradientScaleOnRho(model, ramp, {10, 10}, 0.5, model.layout),
                std::sqrt(5.) / (3. + 10. + 20. + 0.5), 1e-14);
}

TEST(GradientScale, VectorRotationAtConstantMagnitudeIsSeen)
{
    using Strategy = GradientScaleTaggerStrategy<MockMHDModel<2>>;
    MockMHDModel<2> model;
    double constexpr k = 0.3;
    fillRhoV(
        model, [&](auto const& idx) { return std::cos(k * idx[0]); },
        [&](auto const& idx) { return std::sin(k * idx[0]); }, [](auto const&) { return 0.; });
    EXPECT_NEAR((vectorIndicatorAt<Strategy>(model.state.rhoV, std::array{10u, 10u},
                                             GradientScaleParameters{0.1}, model.layout)),
                std::sin(k) / (1. + 0.1), 1e-14);
}

TEST(GradientScale, AnisotropicMeshDoesNotChangeValue)
{
    using Layout = MockMHDModel<2>::gridlayout_type;
    MockMHDModel<2> model;
    Layout const anisotropic{{0.1, 0.4}, {20u, 20u}, PHARE::core::Point<double, 2>{0., 0.}};
    auto const ramp = [](auto const& idx) { return 3. + idx[0] + 2. * idx[1]; };
    EXPECT_NEAR(gradientScaleOnRho(model, ramp, {10, 10}, 0.5, anisotropic),
                std::sqrt(5.) / (3. + 10. + 20. + 0.5), 1e-14);
}

TEST(GradientScale, ZeroFieldIsZero)
{
    MockMHDModel<2> model;
    EXPECT_EQ(
        gradientScaleOnRho(model, [](auto const&) { return 0.; }, {10, 10}, 0.5, model.layout), 0.);
}

TEST(GradientScale, OneAndThreeDimensions)
{
    {
        MockMHDModel<1> model;
        EXPECT_NEAR(
            gradientScaleOnRho(
                model, [](auto const& idx) { return 2. + 0.5 * idx[0]; }, {10}, 1., model.layout),
            0.5 / (7. + 1.), 1e-14);
    }
    {
        MockMHDModel<3> model;
        EXPECT_NEAR(
            gradientScaleOnRho(
                model,
                [](auto const& idx) { return static_cast<double>(idx[0] + idx[1] + idx[2]); },
                {5, 5, 5}, 1., model.layout),
            std::sqrt(3.) / (15. + 1.), 1e-14);
    }
}

std::string gradientScaleError(PHARE::initializer::PHAREDict const& dict)
{
    try
    {
        GradientScaleTaggerStrategy<MockMHDModel<2>> const strategy{dict};
    }
    catch (std::runtime_error const& e)
    {
        return e.what();
    }
    return {};
}

TEST(GradientScale, MissingRefMagnitudeThrows)
{
    auto const missing = gradientScaleError(quantitiesDict({{"rho", 0.1}}));
    EXPECT_TRUE(contains(missing, "rho") and contains(missing, "ref_magnitude")) << missing;
    for (double const refMagnitude : {0., -1., std::numeric_limits<double>::quiet_NaN(),
                                      std::numeric_limits<double>::infinity()})
    {
        auto const error = gradientScaleError(
            quantitiesDict({{"rho", 0.1}}, {{"rho", {{"ref_magnitude", refMagnitude}}}}));
        EXPECT_TRUE(contains(error, "rho") and contains(error, "ref_magnitude")) << error;
    }
    EXPECT_FALSE(gradientScaleError(PHARE::initializer::PHAREDict{}).empty());
    EXPECT_TRUE(
        gradientScaleError(quantitiesDict({{"rho", 0.1}}, {{"rho", {{"ref_magnitude", 1e-3}}}}))
            .empty());
}

TEST(GradientScale, PerQuantityRefMagnitudes)
{
    using Model    = MockMHDModel<2>;
    using Strategy = GradientScaleTaggerStrategy<Model>;
    Model model;
    auto const ghosts = Model::gridlayout_type::options.field_ghost_width;
    fillMHDStep(model, 1., 2.);
    fillRhoV(
        model, [&](auto const& idx) { return idx[0] < ghosts + 5 ? 0. : 1e-3; },
        [](auto const&) { return 0.; }, [](auto const&) { return 0.; });

    auto const allColumnsBut = [&](std::vector<std::uint32_t> const& columns) {
        std::vector<std::uint32_t> others;
        for (std::uint32_t ix = 0; ix < model.layout.nbrCells()[0]; ++ix)
            if (std::find(columns.begin(), columns.end(), ix) == columns.end())
                others.push_back(ix);
        return others;
    };

    auto const rhoSmall = methodTags<Strategy>(
        model, quantitiesDict({{"rho", 0.2}}, {{"rho", {{"ref_magnitude", 0.1}}}}));
    expectColumns(model, rhoSmall, {9, 10}, allColumnsBut({9, 10}));

    EXPECT_TRUE(
        allEqual(methodTags<Strategy>(
                     model, quantitiesDict({{"rho", 0.2}}, {{"rho", {{"ref_magnitude", 100.}}}})),
                 0));

    auto const rhoVSmall = methodTags<Strategy>(
        model, quantitiesDict({{"rhoV", 0.2}}, {{"rhoV", {{"ref_magnitude", 1e-8}}}}));
    expectColumns(model, rhoVSmall, {4, 5}, allColumnsBut({4, 5}));

    EXPECT_TRUE(
        allEqual(methodTags<Strategy>(
                     model, quantitiesDict({{"rhoV", 0.2}}, {{"rhoV", {{"ref_magnitude", 100.}}}})),
                 0));

    auto const both = methodTags<Strategy>(
        model,
        quantitiesDict({{"rho", 0.2}, {"rhoV", 0.2}},
                       {{"rho", {{"ref_magnitude", 100.}}}, {"rhoV", {{"ref_magnitude", 1e-8}}}}));
    expectColumns(model, both, {4, 5}, allColumnsBut({4, 5}));
}

TEST(TaggerFactory, hybridAndMHDMethods)
{
    using hybrid_model = PHARE::solver::PHARE_Types<PHARE::SimOpts{1ul, 1ul, 2ul}>::Hybrid::Model_t;
    using mhd_model    = PHARE::solver::PHARE_Types<mhdOpts(1)>::MHD::Model_t;

    for (std::string const method : {"default", "lohner-per-direction", "lohner-cross-derivatives"})
    {
        PHARE::initializer::PHAREDict dict;
        dict["hybrid_method"] = method;
        dict["mhd_method"]    = method;
        EXPECT_TRUE(TaggerFactory<hybrid_model>::make(dict) != nullptr) << method;
        EXPECT_TRUE(TaggerFactory<mhd_model>::make(dict) != nullptr) << method;
    }

    {
        auto dict             = quantitiesDict({{"B", 0.1}}, {{"B", {{"ref_magnitude", 0.1}}}});
        dict["hybrid_method"] = std::string{"gradient-scale"};
        dict["mhd_method"]    = std::string{"gradient-scale"};
        EXPECT_TRUE(TaggerFactory<hybrid_model>::make(dict) != nullptr);
        EXPECT_TRUE(TaggerFactory<mhd_model>::make(dict) != nullptr);
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
