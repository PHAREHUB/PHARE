#include "tools/bench/core/bench.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"



template<std::size_t dim, std::size_t interp>
void updater_routine(benchmark::State& state)
{
    using namespace PHARE;
    using namespace PHARE::core;

    constexpr std::uint32_t cells   = 30;
    constexpr std::uint32_t n_parts = 1e7;
    auto static constexpr opts      = SimOpts{dim, interp};

    using PHARE_Types   = PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Electromag_t  = UsableElectromag<dim>;
    using ParticleArray = PHARE_Types::ParticleArray_t;
    using Particle_t    = ParticleArray::value_type;
    using Ions          = UsableIons<ParticleArray, interp>;
    using IonUpdater    = IonUpdaterProxy<IonUpdater2<Ions>>;
    using Boxing_t      = UpdaterSelectionBoxing<GridLayout_t, ParticleArray>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    Ions ions{layout};
    Boxing_t const boxing{layout, {grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts())}};

    auto& patch_particles = ions.populations[0].particles;
    patch_particles.domain_particles.vector()
        = std::vector<Particle_t>(n_parts, bench::particle<dim>());
    bench::disperse(patch_particles.domain_particles, 0, cells - 1);
    std::sort(patch_particles.domain_particles);
    auto particles_copy = patch_particles.domain_particles; // tmp storage between update modes

    initializer::PHAREDict dict;
    dict["pusher"]["name"] = std::string{"modified_boris"};
    IonUpdater ionUpdater_{dict};

    double const current_time = 1.0;
    double const new_time     = 1.005;
    auto const dt             = new_time - current_time;
    while (state.KeepRunningBatch(1)) // while (state.KeepRunning())
    {
        ionUpdater_.updatePopulations(UpdaterMode::domain_only, ions, em, boxing, dt);
        ionUpdater_.updateIons(ions);

        patch_particles.domain_particles = particles_copy;
        auto& pack
            = std::get<4>(ions.getRunTimeResourcesViewList()[0].getCompileTimeResourcesViewList());
        pack.setBuffer(&patch_particles.pack());


        ionUpdater_.updatePopulations(UpdaterMode::all, ions, em, boxing, dt);
        ionUpdater_.updateIons(ions);
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
