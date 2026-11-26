#include "push_bench.hpp"
#include "simulator/simulator_def.hpp"

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    auto static constexpr opts      = PHARE::SimOpts{dim, interp};
    constexpr std::uint32_t cells   = 65;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types   = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using Interpolator  = PHARE::core::Interpolator<dim, interp>;
    using Electromag_t  = PHARE::core::UsableElectromag<dim>;
    using Ions_t        = PHARE_Types::Hybrid::Ions_t;
    using ParticleArray = Ions_t::particle_array_type;
    using Particle_t    = ParticleArray::value_type;
    using Pusher        = PHARE::core::IonUpdater<Ions_t, Electromag_t, GridLayout_t>::Pusher;


    GridLayout_t layout{cells};
    Electromag_t em{layout};

    ParticleArray domainParticles{layout.AMRBox()};
    domainParticles.vector()
        = std::vector<Particle_t>(n_parts, PHARE::core::bench::particle<dim>());
    PHARE::core::bench::disperse(domainParticles, 0, cells - 1, 13337);
    // std::sort(domainParticles);

    ParticleArray tmpDomain{layout.AMRBox()};
    tmpDomain.vector() = domainParticles.vector();

    auto rangeIn  = PHARE::core::makeIndexRange(domainParticles);
    auto rangeOut = PHARE::core::makeIndexRange(tmpDomain);

    Interpolator interpolator;
    Pusher const pusher{/*dt=*/.001, /*charge=*/1, /*mass=*/1, layout};
    auto const no_op = [](auto& particleRange) { return particleRange; };
    while (state.KeepRunningBatch(1))
    {
        pusher.move(rangeIn, rangeOut, em, interpolator, no_op, no_op);
    }
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

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
