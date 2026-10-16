#include "core/numerics/ion_updater/ion_updater.hpp"

#include "simulator/simulator_def.hpp"


#include "tools/bench/core/bench.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"


#include "benchmark/benchmark.h"

using namespace PHARE;


template<std::size_t dim, std::size_t interp>
void updater_routine(benchmark::State& state)
{
    auto static constexpr opts    = PHARE::SimOpts{dim, interp};
    constexpr std::uint32_t cells = 15;
    constexpr std::uint32_t ppc   = 100;

    using PHARE_Types                = core::PHARE_Types<opts>;
    using Hybrid_t                   = PHARE_Types::Hybrid;
    using GridLayout_t               = TestGridLayout<typename Hybrid_t::GridLayout_t>;
    auto constexpr static field_opts = core::TensorFieldOptions<Hybrid_t>{};
    using Electromag_t               = core::UsableElectromag<field_opts>;
    using Ions                       = PHARE::core::UsableIons_t<opts>;
    using ParticleArray              = Ions::particle_array_type;
    using IonUpdater_t
        = std::conditional_t<ParticleArray::layout_mode == core::LayoutMode::AoSMapped,
                             core::IonUpdater<ParticleArray, GridLayout_t>,
                             core::ParallelIonUpdater<ParticleArray, GridLayout_t>>;
    using Boxing_t = IonUpdater_t::Boxing_t;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    Ions ions{layout, "protons"};
    Boxing_t const boxing{layout,
                          {grow(layout.AMRBox(), GridLayout_t::options.particle_ghost_width)}};

    auto& patch_particles = ions.populations[0].particles;
    add_particles(patch_particles.domain_particles, layout.AMRBox(), ppc);

    delta_disperse(patch_particles.domain_particles);
    // std::sort(patch_particles.domain_particles);

    auto particles_copy = patch_particles.domain_particles; // tmp storage between update modes

    IonUpdater_t ionUpdater_{};

    double current_time = 1.0;
    double new_time     = 1.005;
    auto dt             = new_time - current_time;
    while (state.KeepRunningBatch(1)) // while (state.KeepRunning())
    {
        ionUpdater_.updatePopulations(ions, em, boxing, dt, core::UpdaterMode::domain_only);
        ions.update();

        patch_particles.domain_particles = particles_copy;
        auto& pack
            = std::get<4>(ions.getRunTimeResourcesViewList()[0].getCompileTimeResourcesViewList());
        pack.setBuffer(&patch_particles.pack());

        ionUpdater_.updatePopulations(ions, em, boxing, dt, core::UpdaterMode::all);
        ions.update();
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
