#include "core/utilities/thread_pool.hpp" // defaults to 1 thread, setup during static init!
#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/data/particles/particle_array_appender.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"
#include "tests/amr/resources_manager/test_resources_manager_fixtures.hpp"

#include "benchmark/benchmark.h"

#include <array>
#include <tuple>
#include <string>
#include <cstddef>
#include <unordered_map>

namespace PHARE::core
{
// RUNTIME ENV VAR OVERRIDES
auto static const cells     = get_env_as("PHARE_CELLS", std::uint32_t{20});
auto static const ppc       = get_env_as("PHARE_PPC", std::size_t{100});
auto static const seed      = get_env_as("PHARE_SEED", std::size_t{1067});
auto static const dt        = get_env_as("PHARE_TIMESTEP", double{.001});
auto static const n_threads = get_env_as("PHARE_THREADS", std::size_t{1});
auto static const n_push    = get_env_as("PHARE_PUSHES", std::size_t{5});

bool static const premain = []() {
    ThreadPool::threads_per_pool = n_threads;
    return true;
}();

template<typename ParticleArray, typename GridLayout>
using IonUpdater_t = std::conditional_t<is_tiled(ParticleArray::layout_mode),
                                        ParallelIonUpdater<ParticleArray, GridLayout>,
                                        IonUpdater<ParticleArray, GridLayout>>;


template<std::size_t dim>
auto make_src_ions(auto const& layout)
{
    UsableIons_t<SimOpts{.dimension = dim}> ions{layout, "protons"};
    auto& particles = ions.populations[0].particles;

    auto const box = layout.AMRBox();
    particles.domain_particles.reserve(box.size() * ppc);
    add_particles(particles.domain_particles, box, ppc);
    add_ghost_particles(particles.level_ghost_particles, box, ppc, 1);

    for (auto* p : {&particles.domain_particles, &particles.level_ghost_particles})
    {
        delta_disperse(*p, seed);
        vary_velocity(*p, -6, 6, seed);
    }
    return ions;
}

template<auto opts>
auto make_ions(auto const& layout, auto const& from)
{
    if constexpr (std::is_same_v<UsableIons_t<opts>, std::decay_t<decltype(from)>>)
        return UsableIons_t<opts>{from};
    else
    {
        UsableIons_t<opts> ions{layout, "protons"};

        auto const add_from = [&]<auto type>(auto const& src, auto& dst) {
            reserve(dst, layout.AMRBox(), ppc);
            append_particles<type>(src, dst);
        };

        auto const& src = from.populations[0].particles;
        auto& dst       = ions.populations[0].particles;
        add_from.template operator()<ParticleType::Domain>(src.domain_particles,
                                                           dst.domain_particles);
        add_from.template operator()<ParticleType::LevelGhost>(src.level_ghost_particles,
                                                               dst.level_ghost_particles);
        return ions;
    }
}


template<auto opts, auto updater_mode>
void update(benchmark::State& state)
{
    using Hybrid_t        = PHARE_Types<opts>::Hybrid;
    using GridLayout_t    = Hybrid_t::GridLayout_t;
    using Ions_t          = UsableIons_t<opts>;
    using EM_t            = UsableElectromag<TensorFieldOptions<Hybrid_t>{}>;
    using Patch_t         = aggregate_adapter<test::HybridPatch<Ions_t, EM_t>>;
    using ParticleArray_t = Ions_t::particle_array_type;
    using Updater_t       = IonUpdater_t<ParticleArray_t, GridLayout_t>;
    using Boxing_t        = Updater_t::Boxing_t;

    using SrcGridLayout_t = PHARE_Types<SimOpts{.dimension = opts.dimension}>::Hybrid::GridLayout_t;

    GridLayout_t const layout        = *TestGridLayout<GridLayout_t>{cells};
    SrcGridLayout_t const src_layout = *TestGridLayout<SrcGridLayout_t>{cells};
    auto const src                   = make_src_ions<opts.dimension>(src_layout);
    Boxing_t const boxing{layout, {layout.AMRBox()}};
    std::unordered_map<std::string, Boxing_t> const boxings{{"patch_id", boxing}};

    std::size_t n_particles = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<Patch_t> patches;
        patches.reserve(1);
        auto& patch = patches.emplace_back(layout, make_ions<opts>(layout, src));
        auto& [_layout, ions, em, electromag] = patch.model.state;
        n_particles += ions.populations[0].particles.domain_particles.size();
        state.ResumeTiming();

        if constexpr (not is_tiled(opts.layout_mode))
            Updater_t{}.updatePopulations(ions, electromag, boxing, dt, updater_mode);
        else
        {
            auto const quantities = [&](int i) {
                return std::forward_as_tuple(patches[i].model.state.ions,
                                             patches[i].model.state.electromag);
            };
            auto accessor = test::make_model_level_accessor(patches, quantities);
            Updater_t{}.updatePopulations(accessor, boxings, dt, updater_mode);
        }
    }
    state.SetItemsProcessed(n_particles);
}


template<std::size_t dim, auto layout_mode, auto updater_mode>
void register_update()
{
    auto constexpr opts   = SimOpts{.dimension = dim, .layout_mode = layout_mode};
    using ParticleArray_t = UsableIons_t<opts>::particle_array_type;

    std::string const mode = updater_mode == UpdaterMode::all ? "all" : "domain_only";
    std::string const name = "update/" + std::to_string(dim) + "d/"
                             + std::string{ParticleArray_t::type_id} + "/" + mode;

    benchmark::RegisterBenchmark(name, &update<opts, updater_mode>)
        ->Unit(benchmark::kMillisecond)
        ->Iterations(n_push)
        ->UseRealTime();
}

} // namespace PHARE::core


int main(int argc, char** argv)
{
    using namespace PHARE::core;
    using enum LayoutMode;

    auto constexpr static layouts = std::array{AoSMapped, AoSTS, AoSPCTS, AoSCMTS};

    for_N<3>([](auto d) {
        constexpr std::size_t dim = d() + 1;
        for_N<layouts.size()>([](auto l) {
            register_update<dim, layouts[l()], UpdaterMode::domain_only>();
            register_update<dim, layouts[l()], UpdaterMode::all>();
        });
    });

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
