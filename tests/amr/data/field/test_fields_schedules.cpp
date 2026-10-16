
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/types.hpp"

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"

#include "gtest/gtest.h"

#include <map>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>


namespace PHARE::amr
{

// few particles: the fields are what is measured here, not the particles
static constexpr std::size_t ppc = 1;

// RUNTIME ENV VAR OVERRIDES
auto static const n_repeats = core::get_env_as("PHARE_REPEATS", std::size_t{3});


template<SimOpts opts>
struct TestParam
{
    auto constexpr static dim    = opts.dimension;
    auto constexpr static interp = opts.interp_order;
    auto constexpr static layout = opts.layout_mode;

    using PhareTypes   = core::PHARE_Types<opts>;
    using GridLayout_t = PhareTypes::Hybrid::GridLayout_t;
    using Hierarchy_t  = AfullHybridBasicHierarchy<opts>;
};

template<typename TestParam_>
struct FieldScheduleHierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    using Hierarchy_t         = TestParam::Hierarchy_t;
    auto constexpr static dim = TestParam::dim;

    std::string const configFile
        = "test_fields_schedules_inputs/" + std::to_string(dim) + "d_config.txt";
    Hierarchy_t hierarchy{configFile, ppc};
};


// AoSMapped first: it is the reference the tiled layouts are compared against
// clang-format off
using FieldDatas = testing::Types<
    TestParam<SimOpts{.dimension=3, .layout_mode=LayoutMode::AoSMapped}>
   ,TestParam<SimOpts{.dimension=3, .layout_mode=LayoutMode::AoSPCTS}>
   ,TestParam<SimOpts{.dimension=3, .layout_mode=LayoutMode::AoSCMTS}>
>;
// clang-format on

TYPED_TEST_SUITE(FieldScheduleHierarchyTest, FieldDatas, );


// reference results of the first (AoSMapped) instantiation, per test and dimension
auto& references()
{
    static std::map<std::string, std::vector<double>> refs;
    return refs;
}

void check_against_reference(std::string const& key, std::vector<double> const& values)
{
    auto& refs = references();
    if (!refs.contains(key))
    {
        refs[key] = values;
        return;
    }
    auto const& ref = refs.at(key);
    ASSERT_EQ(ref.size(), values.size()) << key;
    for (std::size_t i = 0; i < ref.size(); ++i)
        EXPECT_NEAR(ref[i], values[i], 1e-10 * std::max(1., std::abs(ref[i])))
            << key << " value " << i;
}


// sets every domain value (tile domains for tiled fields) from its AMR index
template<typename Field_t>
void set_domain(Field_t& field, auto const& layout, auto&& fn)
{
    if constexpr (core::is_field_tile_set_v<Field_t>)
    {
        for (auto& tile : field())
        {
            auto const tile_gb = tile.ghost_box();
            for (auto const& bix : tile.field_box())
                tile()((bix - tile_gb.lower).as_unsigned()) = fn(bix);
        }
    }
    else
    {
        auto const ghost_box = layout.AMRGhostBoxFor(field);
        for (auto const& bix : layout.AMRBoxFor(field))
            field((bix - ghost_box.lower).as_unsigned()) = fn(bix);
    }
}

// weighted sums over either the ghost layer or the domain of the (reduced) field
// the index weight makes the sum sensitive to values landing in the wrong cell
enum class Region { Ghost, Domain };
template<Region region>
double checksum(auto const& field, auto const& layout, bool& has_nan)
{
    auto constexpr static dim = std::decay_t<decltype(layout)>::dimension;
    auto const& grid          = core::reduce_single(field);
    auto const ghost_box  = layout.AMRGhostBoxFor(field);
    auto const domain_box = layout.AMRBoxFor(field);

    double sum = 0;
    for (auto const& bix : ghost_box)
    {
        if ((region == Region::Ghost) == isIn(bix, domain_box))
            continue;
        auto const lix = (bix - ghost_box.lower).as_unsigned();
        auto const v   = grid(lix);
        has_nan |= std::isnan(v);
        sum += v * (1 + (lix[0] * 7 + lix[dim - 1] * 13) % 31);
    }
    return sum;
}

template<typename Fn>
double time_ms(Fn&& fn)
{
    auto const start = std::chrono::steady_clock::now();
    fn();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
        .count();
}



TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_border_sum_schedules)
{
    auto constexpr static dim = TypeParam::dim;
    using GridLayout_t        = TypeParam::GridLayout_t;
    using Interpolating_t     = core::Interpolating<dim, TypeParam::interp, /*atomic*/ false>;

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    Interpolating_t interpolate;
    for (auto& patch : *lvl0)
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, ions);
        core::resetMoments(ions);
        core::depositParticles(ions, layout, interpolate, core::DomainDeposit{});
    }

    this->hierarchy.messenger->fillDensityBorders(ions, *lvl0, 0);
    this->hierarchy.messenger->fillFluxBorders(ions, *lvl0, 0);

    // particle positions/velocities are not seeded, but particle weights only depend on the
    // density profile, so the summed domain density after border sums must match AoSMapped
    double total  = 0;
    bool has_nan  = false;
    for (auto& patch : *lvl0)
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            auto const& grid      = core::reduce_single(pop.particleDensity());
            auto const ghost_box  = layout.AMRGhostBoxFor(pop.particleDensity());
            auto const domain_box = layout.AMRBoxFor(pop.particleDensity());
            for (auto const& bix : domain_box)
            {
                auto const v = grid((bix - ghost_box.lower).as_unsigned());
                has_nan |= std::isnan(v);
                total += v;
            }
        }
    }

    EXPECT_FALSE(has_nan);
    check_against_reference("border_sum_" + std::to_string(dim), {total});
}



TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_field_refine_schedules)
{
    auto constexpr static dim = TypeParam::dim;
    using GridLayout_t        = TypeParam::GridLayout_t;

    auto& hier      = *this->hierarchy.basicHierarchy->hierarchy();
    auto& rm        = *this->hierarchy.resourcesManagerHybrid;
    auto& model     = *this->hierarchy.hybridModel;
    auto& messenger = *this->hierarchy.messenger;
    auto& em        = model.state.electromag;

    ASSERT_EQ(hier.getNumberOfLevels(), 2);
    auto lvl0 = hier.getPatchLevel(0);
    auto lvl1 = hier.getPatchLevel(1);

    // deterministic, distinct per component and level
    auto const value_fn = [](int const lvl, int const comp) {
        return [=](auto const& bix) {
            double v = 1 + lvl + comp * .5;
            for (std::size_t i = 0; i < dim; ++i)
                v += (i + 1) * .01 * bix[i];
            return v;
        };
    };

    for (int ilvl = 0; ilvl < 2; ++ilvl)
        for (auto& patch : rm.enumerate(*hier.getPatchLevel(ilvl), model))
        {
            auto const layout = layoutFromPatch<GridLayout_t>(*patch);
            for (int c = 0; c < 3; ++c)
            {
                set_domain(em.E[c], layout, value_fn(ilvl, c));
                set_domain(em.B[c], layout, value_fn(ilvl, c + 3));
            }
        }

    // L0 first, its ghosts are the coarse data for L1's level border
    messenger.fillElectricGhosts(em.E, *lvl0, 0);
    messenger.fillMagneticGhosts(em.B, *lvl0, 0);

    double e_ms = 0, b_ms = 0;
    for (std::size_t i = 0; i < n_repeats; ++i)
    {
        e_ms += time_ms([&]() { messenger.fillElectricGhosts(em.E, *lvl1, 0); });
        b_ms += time_ms([&]() { messenger.fillMagneticGhosts(em.B, *lvl1, 0); });
    }

    std::size_t n_cells = 0;
    for (auto& patch : *lvl1)
        n_cells += patch->getBox().size();

    PHARE_LOG_LINE_SS(core::enum_name(TypeParam::layout)
                      << " L1 patches " << lvl1->getLocalNumberOfPatches() << " cells " << n_cells
                      << " fillElectricGhosts " << e_ms / n_repeats << " ms"
                      << " fillMagneticGhosts " << b_ms / n_repeats << " ms");

    std::vector<double> sums;
    bool has_nan = false;
    for (auto& patch : rm.enumerate(*lvl1, model))
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        for (auto const* vf : {&em.E, &em.B})
            for (auto const& field : *vf)
                sums.push_back(checksum<Region::Ghost>(field, layout, has_nan));
    }

    EXPECT_FALSE(has_nan) << "NaN left in L1 ghost cells";
    check_against_reference("field_refine_" + std::to_string(dim), sums);
}



} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
