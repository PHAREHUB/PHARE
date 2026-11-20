#include "push_bench.hpp"
#include "simulator/simulator_def.hpp"
#include "core/numerics/ion_updater/reference/ion_updater.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "benchmark/benchmark.h"


template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    using namespace PHARE::core;

    auto static constexpr opts      = PHARE::SimOpts{dim, interp};
    constexpr std::uint32_t cells   = 65;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types   = PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using Electromag_t  = UsableElectromag<dim>;
    using Interpolator  = PHARE::core::Interpolator<dim, interp>;
    using ParticleArray = PHARE_Types::Hybrid::ParticleArray_t;
    using Particle_t    = ParticleArray::value_type;
    using Ions_t        = PHARE_Types::Hybrid::Ions_t;
    using Pusher        = IonUpdater<Ions_t, Electromag_t, GridLayout_t>::Pusher;
    using Boxing_t      = UpdaterSelectionBoxing<GridLayout_t, ParticleArray>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    UsableIons_t<ParticleArray, interp> ions{layout};

    auto& domainParticles    = ions[0].domainParticles();
    domainParticles.vector() = std::vector<Particle_t>(n_parts, bench::particle<dim>());
    bench::disperse(domainParticles, 0, cells - 1, 13337);

    Boxing_t const boxing{layout,
                          {grow(layout.AMRBox(), GridLayout_t::options.particle_ghost_width)}};
    Interpolator interpolator;
    Pusher const pusher{/*dt=*/.001, /*charge=*/1, /*mass=*/1, layout};
    while (state.KeepRunningBatch(1))
        pusher.template move_and_deposit</*copy_particle=*/true>(ions[0], em, boxing, interpolator);
}


template<std::size_t dim, std::size_t interp>
void push_reference(benchmark::State& state)
{
    using namespace PHARE::core;

    auto static constexpr opts      = PHARE::SimOpts{dim, interp};
    constexpr std::uint32_t cells   = 65;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types   = PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using Electromag_t  = UsableElectromag<dim>;
    using Interpolator  = PHARE::core::Interpolator<dim, interp>;
    using Ions_t        = PHARE_Types::Hybrid::Ions_t;
    using ParticleArray = Ions_t::particle_array_type;
    using Particle_t    = ParticleArray::value_type;
    using Pusher        = reference::IonUpdater<Ions_t, Electromag_t, GridLayout_t>::Pusher;

    GridLayout_t layout{cells};
    Electromag_t em{layout};

    ParticleArray domainParticles{layout.AMRBox()};
    domainParticles.vector() = std::vector<Particle_t>(n_parts, bench::particle<dim>());
    bench::disperse(domainParticles, 0, cells - 1, 13337);

    ParticleArray tmpDomain{layout.AMRBox()};
    tmpDomain.vector() = domainParticles.vector();

    auto rangeIn  = makeIndexRange(domainParticles);
    auto rangeOut = makeIndexRange(tmpDomain);

    Interpolator interpolator;
    Pusher const pusher{/*dt=*/.001, /*charge=*/1, /*mass=*/1, layout};
    auto const no_op = [](auto& particleRange) { return particleRange; };
    while (state.KeepRunningBatch(1))
        pusher.move(rangeIn, rangeOut, em, interpolator, no_op, no_op);
}


BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push_reference, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push_reference, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push_reference, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push_reference, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push_reference, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push_reference, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push_reference, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push_reference, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push_reference, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
