
#include <stdexcept>
#include "core/logger.hpp"

#include "core/utilities/thread_pool.hpp" // defaults to 1 thread, setup during static init!
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp" // IWYU pragma: keep
#include "core/data/particles/particle_array_appender.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/grid/test_grid_fixtures.hpp"
#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"
#include "tests/amr/resources_manager/test_resources_manager_fixtures.hpp"

#include "gtest/gtest.h"

#include <tuple>
#include <cstddef>

namespace PHARE::core
{
// RUNTIME ENV VAR OVERRIDES
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{30});
auto static const ppc       = get_env_as("PHARE_PPC", std::size_t{100});
auto static const seed      = get_env_as("PHARE_SEED", std::size_t{1067});
auto static const n_patches = get_env_as("PHARE_PATCHES", std::size_t{1});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const n_threads = get_env_as("PHARE_THREADS", std::size_t{1});
auto static const n_push    = get_env_as("PHARE_PUSHES", std::size_t{1});

bool static const premain = []() {
    PHARE_WITH_PHLOP({
        PHARE_LOG_LINE_STR("cells      : " << cells);
        PHARE_LOG_LINE_STR("ppc        : " << ppc);
        PHARE_LOG_LINE_STR("particles  : " << std::pow(cells, 3) * ppc);
        PHARE_LOG_LINE_STR("n_patches  : " << n_patches);
        PHARE_LOG_LINE_STR("seed       : " << seed);
        PHARE_LOG_LINE_SS("particle MB ≈ " << n_patches * std::pow(cells, 3) * ppc * 76 / 1e6);
    })
    phlop::scope_timer().file_name(".phare_times.0.txt").init();
    ThreadPool::threads_per_pool = n_threads;
    return true;
}();

auto& pool = *ThreadPool::INSTANCE().thread_pools[0];

template<typename Patches>
void ref_update(UpdaterMode mode, Patches& patches)
{
    using GridLayout_t = Patches::value_type::GridLayout_t;

    for (std::size_t i = 0; i < n_push; ++i)
        for (auto& patch : patches)
        {
            auto& [layout, ions, em, electromag] = patch.model.state;
            using IonUpdater_t = IonUpdater<std::decay_t<decltype(ions)>, std::decay_t<decltype(electromag)>,
                                            GridLayout_t>;
            using Boxing_t
                = std::decay_t<decltype(*selection_boxing_impl<IonUpdater_t, GridLayout_t>())>;

            Boxing_t const boxing{layout, {layout.AMRBox()}};
            IonUpdater_t{}.updatePopulations(ions, electromag, boxing, dt, mode);
            ions.computeChargeDensity();
            ions.computeBulkVelocity();
        }
}

template<typename Patches>
void cmp_update(UpdaterMode mode, Patches& patches)
{
    using GridLayout_t = Patches::value_type::GridLayout_t;
    using Boxing_t     = PHARE::core::UpdaterSelectionBoxing<GridLayout_t>;

    auto const& layout = patches[0].layout;
    std::unordered_map<std::string, Boxing_t> boxings;
    boxings.try_emplace("patch_id", Boxing_t{layout, {layout.AMRBox()}});

    auto const quantities = [&](int i) {
        return std::forward_as_tuple(patches[0].model.state.ions,
                                     patches[i].model.state.electromag);
    };
    auto accessor = test::make_model_level_accessor(patches, quantities);

    using Ions_t       = std::decay_t<decltype(patches[0].model.state.ions)>;
    using Electromag_t = std::decay_t<decltype(patches[0].model.state.electromag)>;

    try
    {
        for (std::size_t i = 0; i < n_push; ++i)
        {
            IonUpdater<Ions_t, Electromag_t, GridLayout_t, decltype(accessor)>{accessor}
                .updatePopulations(boxings, dt, mode);

            for (auto& patch : patches)
            {
                auto& [layout, ions, em, electromag] = patch.model.state;
                ions.computeChargeDensity();
                ions.computeBulkVelocity();
            }
        }
    }
    catch (std::exception const& e)
    {
        PHARE_LOG_LINE_SS(e.what());
    }
}


template<auto opts>
auto make_ref_ions(auto const& layout)
{
    UsableIons_t<opts> ions{layout, "protons"};

    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto const disperse = [&](auto& particles) {
        delta_disperse(particles.domain_particles, seed);
        vary_velocity(particles.domain_particles, -6, 6, seed);
    };

    auto const particle_box = layout.AMRBox();
    auto add_particles      = [&](auto& particles) {
        particles.domain_particles.reserve(particle_box.size() * ppc);
        core::add_particles(particles.domain_particles, particle_box, ppc);
    };

    add_particles(ions.populations[0].particles);
    disperse(ions.populations[0].particles);


    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), particle_box.size() * ppc);

    return ions;
}


template<auto opts>
auto from_ions(auto const& layout, auto const& from)
{
    UsableIons_t<opts> ions{layout, "protons"};
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), 0ull);

    auto _add_particles_from = [&]<auto type>(auto& src, auto& dst) {
        reserve(dst, layout.AMRBox(), ppc);
        append_particles<type>(src, dst);
    };

    _add_particles_from.template operator()<ParticleType::Domain>(
        from.populations[0].particles.domain_particles,
        ions.populations[0].particles.domain_particles);

    auto const particle_box = layout.AMRBox();
    EXPECT_EQ(ions.populations[0].particles.domain_particles.size(), particle_box.size() * ppc);
    return ions;
}


template<auto type, typename GridLayout_t, typename P0, typename P1>
void check_particles(GridLayout_t const& layout, P0& ref, P1& cmp_, double const atol,
                     std::string const& name)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    using CPU_ref = ParticleArray<ParticleArrayOptions{P0::dimension, LayoutMode::AoS,
                                                       StorageMode::VECTOR, AllocatorMode::CPU}>;

    auto const box = layout.AMRBox();
    // level/patch ghost particles live outside the domain box entirely - the sorter's
    // flat-index cell key isn't defined out there, so sorting with the plain domain box
    // can scatter/misgroup ghost-region particles differently between ref and cmp even
    // when the underlying data matches. Grow the box for non-domain comparisons so the
    // flattener covers every cell these particles can actually occupy.
    auto const sort_box = type == ParticleType::Domain
                              ? box
                              : grow(box, GridLayout_t::options.particle_ghost_width);

    auto cmp = convert_particles<CPU_ref>(cmp_, layout);
    sort_particles(cmp, sort_box);
    sort_particles(ref, sort_box);

    EXPECT_EQ(ref.size(), cmp.size());

    auto const report = compare_particles(ref, cmp, atol);
    if (report)
    {
        PHARE_LOG_LINE_STR("Comparing Particle Arrays OK: " << P0::id() << " vs " << P1::id());
    }
    else
    {
        PHARE_LOG_LINE_STR("Comparing Particle Arrays FAIL: " << P0::id() << " vs " << P1::id());
    }
    PHARE_LOG_LINE_STR("results: " << report.why());
    if (ref.size())
        PHARE_LOG_LINE_STR("eg: " << ref[0]);

    EXPECT_TRUE(report);
}

template<typename GridLayout_t, typename R, typename C>
void compare(GridLayout_t const& layout, R& ref, C& cmp)
{
    using ParticleArray_t = C::ParticleArray_t;

    using enum LayoutMode;
    using enum AllocatorMode;

    double diff = 1e-15;
    if constexpr (is_tiled(ParticleArray_t::layout_mode))
        diff *= 1e1; // p2m op order diff

    check_particles<ParticleType::Domain>(layout, ref.populations[0].particles.domain_particles,
                                          cmp.populations[0].particles.domain_particles, diff,
                                          "domain");

    check_particles<ParticleType::PatchGhost>(
        layout, ref.populations[0].particles.patch_ghost_particles,
        cmp.populations[0].particles.patch_ghost_particles, diff, "patchghost");

    {
        // zero_ghost_layer(layout, ref.rhoC, cmp.rhoC);
        auto const rhoCport = compare_reduced_fields(ref.rhoC, cmp.rhoC, diff);
        PHARE_LOG_LINE_STR("results: ions rhoC" << rhoCport.why());
        EXPECT_TRUE(rhoCport);

        // zero_ghost_layer(layout, ref.rhoM, cmp.rhoM);
        auto const rhoMport = compare_reduced_fields(ref.rhoM, cmp.rhoM, diff);
        PHARE_LOG_LINE_STR("results: ions rhoM" << rhoMport.why());
        EXPECT_TRUE(rhoMport);
    }

    auto const freport
        = compare_reduced_tensor_fields(ref.populations[0].F, cmp.populations[0].F, diff);
    PHARE_LOG_LINE_STR("results: " << freport.why());
    EXPECT_TRUE(freport);

    auto const rhoCport
        = compare_reduced_fields(ref.populations[0].rhoC, cmp.populations[0].rhoC, diff);
    PHARE_LOG_LINE_STR("results: " << rhoCport.why());
    EXPECT_TRUE(rhoCport);

    auto const rhoPport
        = compare_reduced_fields(ref.populations[0].rhoP, cmp.populations[0].rhoP, diff);
    PHARE_LOG_LINE_STR("results: " << rhoPport.why());
    EXPECT_TRUE(rhoPport);
}


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode, auto _updater_mode>
struct TestParam
{
    static_assert(std::is_same_v<decltype(_layout_mode), LayoutMode>);
    static_assert(std::is_same_v<decltype(_alloc_mode), AllocatorMode>);
    auto constexpr static dim          = _dim;
    auto constexpr static layout_mode  = _layout_mode;
    auto constexpr static alloc_mode   = _alloc_mode;
    auto constexpr static updater_mode = _updater_mode;
    auto constexpr static opts
        = SimOpts{.dimension = _dim, .layout_mode = _layout_mode, .alloc_mode = _alloc_mode};
};



template<typename Param>
struct MultiPatchIonUpdaterTest : public ::testing::Test
{
    auto constexpr static dim          = Param::dim;
    auto constexpr static layout_mode  = Param::layout_mode;
    auto constexpr static alloc_mode   = Param::alloc_mode;
    auto constexpr static updater_mode = Param::updater_mode;

    auto constexpr static ref_opts = SimOpts{.dimension = dim};
    auto constexpr static cmp_opts = Param::opts;

    using RefCoreTypes = PHARE_Types<ref_opts>::Hybrid;
    using CmpCoreTypes = PHARE_Types<cmp_opts>::Hybrid;

    using RefGridLayout_t = RefCoreTypes::GridLayout_t;
    using CmpGridLayout_t = CmpCoreTypes::GridLayout_t;

    using TestRefGridLayout_t = TestGridLayout<RefGridLayout_t>;
    using TestCmpGridLayout_t = TestGridLayout<CmpGridLayout_t>;

    using RefParticleArray_t = AoSMappedParticleArray<dim>;
    using CmpParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    auto constexpr static ref_field_opts = TensorFieldOptions<RefCoreTypes>{};
    auto constexpr static cmp_field_opts = TensorFieldOptions<CmpCoreTypes>{};

    using RefIons_t = UsableIons_t<ref_opts>;
    using RefEM_t   = UsableElectromag<ref_field_opts>;
    using CmpIons_t = UsableIons_t<cmp_opts>;
    using CmpEM_t   = UsableElectromag<cmp_field_opts>;

    RefGridLayout_t const ref_layout = *TestRefGridLayout_t{cells};
    CmpGridLayout_t const cmp_layout = *TestCmpGridLayout_t{cells};

    MultiPatchIonUpdaterTest() {}

    void run()
    {
        cmp_patches.reserve(n_patches);
        auto& ref = ref_patches.emplace_back(ref_layout, make_ref_ions<ref_opts>(ref_layout));

        for (std::size_t i = 0; i < n_patches; i++)
            cmp_patches.emplace_back(cmp_layout,
                                     from_ions<cmp_opts>(cmp_layout, ref.model.state.ions));

        ref_update(updater_mode, ref_patches);
        cmp_update(updater_mode, cmp_patches);

        for (auto& cmp : cmp_patches)
            compare(cmp_layout, ref_patches[0].model.state.ions, cmp.model.state.ions);
    }



    using RefPatch = test::HybridPatch<RefIons_t, RefEM_t>;
    using CmpPatch = test::HybridPatch<CmpIons_t, CmpEM_t>;

    std::vector<aggregate_adapter<RefPatch>> ref_patches{};
    std::vector<aggregate_adapter<CmpPatch>> cmp_patches{};
};

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !
    TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU, UpdaterMode::domain_only>
   ,TestParam<1, LayoutMode::AoSPCTS, AllocatorMode::CPU, UpdaterMode::all>
   ,TestParam<2, LayoutMode::AoSPCTS, AllocatorMode::CPU, UpdaterMode::domain_only>
   ,TestParam<2, LayoutMode::AoSPCTS, AllocatorMode::CPU, UpdaterMode::all>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU, UpdaterMode::domain_only>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU, UpdaterMode::all>

>;
// clang-format on

TYPED_TEST_SUITE(MultiPatchIonUpdaterTest, Permutations_t, );

TYPED_TEST(MultiPatchIonUpdaterTest, test_updater)
{
    this->run();
}

template<auto opts> // used by gtest
void PrintTo(ParticleArray<opts> const& arr, std::ostream* os)
{
    *os << arr;
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    try
    {
        auto r = RUN_ALL_TESTS();
        PHARE_WITH_PHLOP(phlop::scope_timer().shutdown());
        return r;
    }
    catch (std::runtime_error const& e)
    {
        PHARE_LOG_LINE_SS(e.what());
    }
    return 1;
}
