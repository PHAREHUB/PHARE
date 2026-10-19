

#include "phare_solver.hpp"
#include "simulator/simulator_def.hpp"

#include "amr/tagging/tagger_factory.hpp"
#include "amr/tagging/tagging_quantities.hpp"
#include "amr/tagging/default_tagger_strategy.hpp"
#include "amr/tagging/lohner_tagger_strategy.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"
#include "tests/core/data/mhd_state/test_mhd_state_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <limits>
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

template<std::size_t dim>
double legacyProfile(std::size_t component, std::array<std::uint32_t, dim> const& idx)
{
    if (dim == 1 and component == 0)
        return 1.0;
    double wave = 0.;
    for (std::size_t d = 1; d < dim; ++d)
        wave += std::sin(0.9 * idx[d] + 0.4 * d + component);
    return (1. + component) * std::tanh((idx[0] - 11.0) / 2.0) + 0.05 * wave;
}


template<typename GridLayoutT, typename Field>
void legacyDefaultTags(Field const& Bx, Field const& By, Field const& Bz, GridLayoutT const& layout,
                       double threshold_, int* tags)
{
    using gridlayout_type           = GridLayoutT;
    static auto constexpr dimension = GridLayoutT::dimension;

    auto const& start_x
        = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);
    auto const& end_x = layout.nbrCells()[0] - 1;

    bool constexpr c_ordering = false;
    auto tagsv = PHARE::core::NdArrayView<dimension, int, c_ordering>(tags, layout.nbrCells());

    if constexpr (dimension == 1)
    {
        auto constexpr doLastCell = gridlayout_type::options.field_ghost_width > 2;
        std::size_t oneOrZero     = doLastCell ? 1 : 0;
        for (auto iCell = 0u, ix = start_x; iCell < end_x + oneOrZero; ++ix, ++iCell)
        {
            auto Byavg     = 0.2 * (By(ix - 2) + By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2));
            auto Bzavg     = 0.2 * (Bz(ix - 2) + Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2));
            auto Byavgp1   = 0.2 * (By(ix - 1) + By(ix) + By(ix + 1) + By(ix + 2) + By(ix + 3));
            auto Bzavgp1   = 0.2 * (Bz(ix - 1) + Bz(ix) + Bz(ix + 1) + Bz(ix + 2) + Bz(ix + 3));
            auto criter_by = std::abs(Byavgp1 - Byavg) / (1 + std::abs(Byavg));
            auto criter_bz = std::abs(Bzavgp1 - Bzavg) / (1 + std::abs(Bzavg));
            auto criter_b  = std::sqrt(criter_by * criter_by + criter_bz * criter_bz);
            auto criter    = criter_b;

            if (criter > threshold_)
            {
                tagsv(iCell) = 1;
            }
            else
                tagsv(iCell) = 0;
        }
    }
    if constexpr (dimension == 2)
    {
        auto const& start_y
            = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Y);

        auto const& end_y = layout.nbrCells()[1] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x <= end_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y <= end_y; ++iy, ++iTag_y)
            {
                auto field_diff = [&](auto const& F) //
                {
                    auto const delta_2x = std::abs(F(ix + 2, iy) - F(ix, iy));
                    auto const delta_2y = std::abs(F(ix, iy + 2) - F(ix, iy));
                    auto const delta_x  = std::abs(F(ix + 1, iy) - F(ix, iy));
                    auto const delta_y  = std::abs(F(ix, iy + 1) - F(ix, iy));

                    auto const criter_x = delta_2x / (1 + delta_x);
                    auto const criter_y = delta_2y / (1 + delta_y);

                    return std::make_tuple(criter_x, criter_y);
                };

                auto const& [Bx_x, Bx_y] = field_diff(Bx);
                auto const& [By_x, By_y] = field_diff(By);
                auto const& [Bz_x, Bz_y] = field_diff(Bz);
                auto crit                = std::max({Bx_x, Bx_y, By_x, By_y, Bz_x, Bz_y});

                if (crit > threshold_)
                {
                    tagsv(iTag_x, iTag_y) = 1;
                }
                else
                {
                    tagsv(iTag_x, iTag_y) = 0;
                }
            }
        }
    }
    if constexpr (dimension == 3)
    {
        auto const& start_y
            = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Y);
        auto const& start_z
            = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::Z);

        auto const& end_y = layout.nbrCells()[1] - 1;
        auto const& end_z = layout.nbrCells()[2] - 1;

        for (auto iTag_x = 0u, ix = start_x; iTag_x <= end_x; ++ix, ++iTag_x)
        {
            for (auto iTag_y = 0u, iy = start_y; iTag_y <= end_y; ++iy, ++iTag_y)
            {
                for (auto iTag_z = 0u, iz = start_z; iTag_z <= end_z; ++iz, ++iTag_z)
                {
                    auto const field_diff = [&](auto const& F) {
                        auto const delta_2x = std::abs(F(ix + 2, iy, iz) - F(ix, iy, iz));
                        auto const delta_2y = std::abs(F(ix, iy + 2, iz) - F(ix, iy, iz));
                        auto const delta_2z = std::abs(F(ix, iy, iz + 2) - F(ix, iy, iz));
                        auto const delta_x  = std::abs(F(ix + 1, iy, iz) - F(ix, iy, iz));
                        auto const delta_y  = std::abs(F(ix, iy + 1, iz) - F(ix, iy, iz));
                        auto const delta_z  = std::abs(F(ix, iy, iz + 1) - F(ix, iy, iz));

                        auto const criter_x = delta_2x / (1 + delta_x);
                        auto const criter_y = delta_2y / (1 + delta_y);
                        auto const criter_z = delta_2z / (1 + delta_z);

                        return std::make_tuple(criter_x, criter_y, criter_z);
                    };

                    auto const& [Bx_x, Bx_y, Bx_z] = field_diff(Bx);
                    auto const& [By_x, By_y, By_z] = field_diff(By);
                    auto const& [Bz_x, Bz_y, Bz_z] = field_diff(Bz);
                    auto const crit
                        = std::max({Bx_x, Bx_y, Bx_z, By_x, By_y, By_z, Bz_x, Bz_y, Bz_z});

                    if (crit > threshold_)
                    {
                        tagsv(iTag_x, iTag_y, iTag_z) = 1;
                    }
                    else
                    {
                        tagsv(iTag_x, iTag_y, iTag_z) = 0;
                    }
                }
            }
        }
    }
}


template<typename Model>
void fillLegacyB(Model& model)
{
    using PHARE::core::Component;
    auto constexpr dim = Model::dimension;
    auto& B            = model.state.electromag.B;
    for (std::size_t c = 0; c < 3; ++c)
        fillField(B.getComponent(static_cast<Component>(c)),
                  [&](auto const& idx) { return legacyProfile<dim>(c, idx); });
}

template<typename Model>
std::vector<int> sentinelTags(Model const& model)
{
    return std::vector<int>(PHARE::core::product(model.layout.nbrCells()), -1);
}

template<typename Model>
std::vector<int> legacyTags(Model& model, double threshold)
{
    using PHARE::core::Component;
    auto& B   = model.state.electromag.B;
    auto tags = sentinelTags(model);
    legacyDefaultTags(B.getComponent(Component::X), B.getComponent(Component::Y),
                      B.getComponent(Component::Z), model.layout, threshold, tags.data());
    return tags;
}

template<typename Model>
std::vector<int> strategyTags(Model& model, PHARE::initializer::PHAREDict const& dict)
{
    auto tags = sentinelTags(model);
    DefaultTaggerStrategy<Model>{dict}.tag(model, model.layout, tags.data());
    return tags;
}

PHARE::initializer::PHAREDict quantitiesDict(std::vector<tagging::QuantityThreshold> const& qs)
{
    PHARE::initializer::PHAREDict dict;
    dict["nbr_quantities"] = static_cast<int>(qs.size());
    for (std::size_t i = 0; i < qs.size(); ++i)
    {
        auto const path         = "Q" + std::to_string(i);
        dict[path]["name"]      = qs[i].name;
        dict[path]["threshold"] = qs[i].threshold;
    }
    return dict;
}

void expectSomeTaggedSomeNot(std::vector<int> const& tags)
{
    EXPECT_NE(std::count(tags.begin(), tags.end(), 1), 0);
    EXPECT_NE(std::count(tags.begin(), tags.end(), 0), 0);
}


template<typename Model>
struct DefaultStrategyLegacy : public ::testing::Test
{
    DefaultStrategyLegacy() { fillLegacyB(model); }
    Model model{24};
};

using LegacyModels = testing::Types<MockHybridModel<1, 1>, MockHybridModel<1, 2>,
                                    MockHybridModel<2, 1>, MockHybridModel<3, 1>>;
TYPED_TEST_SUITE(DefaultStrategyLegacy, LegacyModels);

TYPED_TEST(DefaultStrategyLegacy, MatchesLegacyOnB)
{
    PHARE::initializer::PHAREDict dict;
    dict["threshold"]   = 0.1;
    auto const expected = legacyTags(this->model, 0.1);
    expectSomeTaggedSomeNot(expected);
    EXPECT_EQ(strategyTags(this->model, dict), expected);
}

TYPED_TEST(DefaultStrategyLegacy, MatchesLegacyWithExplicitBQuantity)
{
    auto const expected = legacyTags(this->model, 0.1);
    EXPECT_EQ(strategyTags(this->model, quantitiesDict({{"B", 0.1}})), expected);
}

TEST(DefaultStrategy, NaNOrderingMatchesLegacy)
{
    using PHARE::core::Component;
    MockHybridModel<2> model{24};
    fillLegacyB(model);
    model.state.electromag.B.getComponent(Component::Y)(13, 9)
        = std::numeric_limits<double>::quiet_NaN();
    auto const expected = legacyTags(model, 0.1);
    EXPECT_EQ(strategyTags(model, quantitiesDict({{"B", 0.1}})), expected);
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



template<typename Fn>
double lohnerOnDualBy(Fn&& fn, std::uint32_t i, double reltol)
{
    MockHybridModel<1> model;
    fillField(model.state.electromag.B.getComponent(PHARE::core::Component::Y),
              [&](auto const& idx) { return fn(idx[0]); });
    return lohner::indicator<1>(tagging::resolveQuantity(model, "By"),
                                std::array<std::uint32_t, 1>{i}, reltol, 1e-30);
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

TEST(LohnerIndicator, CellCenteredAveragesPrimal)
{
    MockHybridModel<1> model;
    auto& B = model.state.electromag.B;
    fillField(B.getComponent(PHARE::core::Component::X),
              [](auto const& idx) { return static_cast<double>(idx[0]); });
    fillField(B.getComponent(PHARE::core::Component::Y),
              [](auto const& idx) { return static_cast<double>(idx[0]); });
    auto const idx = std::array<std::uint32_t, 1>{7};
    EXPECT_DOUBLE_EQ(lohner::cellCentered<1>(tagging::resolveQuantity(model, "Bx")[0], idx), 7.5);
    EXPECT_DOUBLE_EQ(lohner::cellCentered<1>(tagging::resolveQuantity(model, "By")[0], idx), 7.0);
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
