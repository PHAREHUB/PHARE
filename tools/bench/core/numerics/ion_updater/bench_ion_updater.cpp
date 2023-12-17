#include "tools/bench/core/bench.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

using namespace PHARE;

template<typename ParticleArray_t>
struct PatchState
{
    using box_t = typename ParticleArray_t::box_t;
    PatchState(box_t const& domain_)
        : domain{domain_}
    {
    }
    auto pack() { return &particles_pack; }
    box_t domain;
    ParticleArray_t domain_particles{domain, 0}; // ignore size, see next line
    ParticleArray_t patch_ghost_particles = domain_particles;
    ParticleArray_t level_ghost_particles = domain_particles;
    core::ParticlesPack<ParticleArray_t> particles_pack{&domain_particles, &patch_ghost_particles,
                                                        &level_ghost_particles,
                                                        /*levelGhostParticlesOld=*/nullptr,
                                                        /*levelGhostParticlesNew=*/nullptr};
};

template<std::size_t dim, std::size_t interp>
void updater_routine(benchmark::State& state)
{
    constexpr std::uint32_t cells   = 30;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types   = core::PHARE_Types<dim, interp>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Electromag_t  = core::bench::Electromag<GridLayout_t>;
    using Ions          = typename PHARE_Types::Ions_t;
    using IonPop        = typename Ions::value_type;
    using ParticleArray = typename PHARE_Types::ParticleArray_t;
    using Particle_t    = typename ParticleArray::value_type;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    core::bench::Flux<GridLayout_t> flux{layout, "pop0_flux"};
    core::bench::BulkV<GridLayout_t> bulkV{layout};
    auto rho = core::bench::rho(layout);

    PatchState<ParticleArray> patch_particles{layout.AMRBox()};
    patch_particles.domain_particles.vector()
        = std::vector<Particle_t>(n_parts, core::bench::particle<dim>());
    core::bench::disperse(patch_particles.domain_particles, 0, cells - 1);
    std::sort(patch_particles.domain_particles);
    auto particles_copy = patch_particles.domain_particles; // tmp storage between update modes

    std::string popName = "pop0";
    auto ions           = core::bench::single_pop_ions_from<Ions>(popName, rho, bulkV, flux,
                                                        patch_particles.pack());

    initializer::PHAREDict dict;
    dict["pusher"]["name"] = std::string{"modified_boris"};
    core::IonUpdater<Ions, Electromag_t, GridLayout_t> ionUpdater_{dict};

    double current_time = 1.0;
    double new_time     = 1.005;
    auto dt             = new_time - current_time;
    while (state.KeepRunningBatch(1)) // while (state.KeepRunning())
    {
        ionUpdater_.updatePopulations(ions, em, layout, dt, core::UpdaterMode::domain_only);
        ionUpdater_.updateIons(ions, layout);

        patch_particles.domain_particles = particles_copy;
        ions.getRunTimeResourcesUserList()[0].setBuffer(popName, patch_particles.pack());

        ionUpdater_.updatePopulations(ions, em, layout, dt, core::UpdaterMode::all);
        ionUpdater_.updateIons(ions, layout);
    }
}

BENCHMARK_TEMPLATE(updater_routine, 1, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(updater_routine, 1, 2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(updater_routine, 1, 3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(updater_routine, 2, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(updater_routine, 2, 2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(updater_routine, 2, 3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(updater_routine, 3, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(updater_routine, 3, 2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(updater_routine, 3, 3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
