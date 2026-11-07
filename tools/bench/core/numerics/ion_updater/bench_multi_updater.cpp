// push + deposit of one patch, AoSMapped (IonUpdater) vs a tiled layout (ParallelIonUpdater)
//  same particles in both, see tests/core/numerics/ion_updater/test_multi_updater.cpp
//  env: PHARE_CELLS (default 30), PHARE_PPC (default 33)

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

#include "benchmark/benchmark.h"

#include <cstddef>
#include <unordered_map>

namespace PHARE::core
{
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{30});
auto static const ppc   = get_env_as("PHARE_PPC", std::size_t{33});
auto static const seed  = std::size_t{1067};
auto static const dt    = double{.001};

template<auto opts>
auto make_ions(auto const& layout)
{
    UsableIons_t<opts> ions{layout, "protons"};
    auto& particles         = ions.populations[0].particles;
    auto const particle_box = layout.AMRBox();
    reserve(particles.domain_particles, particle_box, ppc);
    core::add_particles(particles.domain_particles, particle_box, ppc);
    add_ghost_particles(particles.level_ghost_particles, particle_box, ppc, 1);
    delta_disperse(particles.domain_particles, seed);
    vary_velocity(particles.domain_particles, -6, 6, seed);
    delta_disperse(particles.level_ghost_particles, seed);
    vary_velocity(particles.level_ghost_particles, -6, 6, seed);
    return ions;
}

template<auto opts>
auto from_ions(auto const& layout, auto const& from)
{
    UsableIons_t<opts> ions{layout, "protons"};
    auto& dst = ions.populations[0].particles;
    auto& src = from.populations[0].particles;
    reserve(dst.domain_particles, layout.AMRBox(), ppc);
    append_particles<ParticleType::Domain>(src.domain_particles, dst.domain_particles);
    append_particles<ParticleType::LevelGhost>(src.level_ghost_particles,
                                               dst.level_ghost_particles);
    return ions;
}

template<std::size_t dim>
struct Setup
{
    auto constexpr static ref_opts = SimOpts{.dimension = dim};
    auto constexpr static cmp_opts = SimOpts{.dimension = dim, .layout_mode = LayoutMode::AoSCMTS};

    using RefTypes        = PHARE_Types<ref_opts>::Hybrid;
    using CmpTypes        = PHARE_Types<cmp_opts>::Hybrid;
    using RefGridLayout_t = RefTypes::GridLayout_t;
    using CmpGridLayout_t = CmpTypes::GridLayout_t;
    using RefPatch        = test::HybridPatch<UsableIons_t<ref_opts>,
                                              UsableElectromag<TensorFieldOptions<RefTypes>{}>>;
    using CmpPatch        = test::HybridPatch<UsableIons_t<cmp_opts>,
                                              UsableElectromag<TensorFieldOptions<CmpTypes>{}>>;

    static Setup& INSTANCE()
    {
        static Setup self;
        return self;
    }

    RefGridLayout_t const ref_layout = *TestGridLayout<RefGridLayout_t>{cells};
    CmpGridLayout_t const cmp_layout = *TestGridLayout<CmpGridLayout_t>{cells};
    std::vector<aggregate_adapter<RefPatch>> ref{};
    std::vector<aggregate_adapter<CmpPatch>> cmp{};

    Setup()
    {
        auto& r = ref.emplace_back(ref_layout, make_ions<ref_opts>(ref_layout));
        cmp.emplace_back(cmp_layout, from_ions<cmp_opts>(cmp_layout, r.model.state.ions));
    }
};


template<std::size_t dim, auto mode>
void ref_update(benchmark::State& state)
{
    auto& patches = Setup<dim>::INSTANCE().ref;
    using GridLayout_t = Setup<dim>::RefGridLayout_t;
    for (auto _ : state)
        for (auto& patch : patches)
        {
            auto& [layout, ions, em, electromag] = patch.model.state;
            using ParticleArray_t = std::decay_t<decltype(ions)>::particle_array_type;
            using Updater_t       = IonUpdaterFor_t<ParticleArray_t, GridLayout_t>;
            typename Updater_t::Boxing_t const boxing{layout, {layout.AMRBox()}};
            Updater_t{}.updatePopulations(ions, electromag, boxing, dt, mode);
        }
}

template<std::size_t dim, auto mode>
void cmp_update(benchmark::State& state)
{
    auto& patches      = Setup<dim>::INSTANCE().cmp;
    using GridLayout_t = Setup<dim>::CmpGridLayout_t;
    using Boxing_t     = UpdaterSelectionBoxing<GridLayout_t>;
    auto const& layout = patches[0].layout;
    std::unordered_map<std::string, Boxing_t> boxings;
    boxings.try_emplace("patch_id", Boxing_t{layout, {layout.AMRBox()}});
    auto const quantities = [&](int i) {
        return std::forward_as_tuple(patches[i].model.state.ions,
                                     patches[i].model.state.electromag);
    };
    auto accessor         = test::make_model_level_accessor(patches, quantities);
    using ParticleArray_t = std::decay_t<decltype(patches[0].model.state.ions)>::particle_array_type;
    for (auto _ : state)
        IonUpdaterFor_t<ParticleArray_t, GridLayout_t>{}.updatePopulations(accessor, boxings, dt,
                                                                           mode);
}

BENCHMARK_TEMPLATE(ref_update, 3, UpdaterMode::domain_only)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(cmp_update, 3, UpdaterMode::domain_only)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(ref_update, 3, UpdaterMode::all)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(cmp_update, 3, UpdaterMode::all)->Unit(benchmark::kMillisecond);

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
