
#include "benchmark/benchmark.h"

#include "pusher_bench.h"

using namespace PHARE::core::bench;

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    constexpr std::uint32_t cells = 65;
    constexpr std::uint32_t parts = 1e8;

    using PHARE_Types       = PHARE::core::PHARE_Types<dim, interp>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Ions_t            = typename PHARE_Types::Ions_t;
    using Electromag_t      = typename PHARE_Types::Electromag_t;
    using GridLayout_t      = typename PHARE_Types::GridLayout_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using PartIterator      = typename ParticleArray::iterator;


    using BorisPusher_t = PHARE::core::BorisPusher<dim, PartIterator, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;

    Interpolator interpolator;
    ParticleArray domainParticles{parts, particle<dim>(/*icell =*/34)};

    auto range    = PHARE::core::makeRange(domainParticles);
    auto meshSize = PHARE::core::ConstArray<double, dim>(1.0 / cells);
    auto nCells   = PHARE::core::ConstArray<std::uint32_t, dim>(cells);
    auto origin   = PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)};

    GridLayout_t layout{meshSize, nCells, origin};

    PHARE::core::bench::Electromag<GridLayout_t, VecField<dim>> electromag{layout};

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    while (state.KeepRunning())
    {
        pusher.move(
            /*ParticleRange const&*/ range,
            /*ParticleRange&*/ range,
            /*Electromag const&*/ electromag,
            /*double mass*/ 1,
            /*Interpolator&*/ interpolator,
            /*ParticleSelector const&*/ [](auto const& /*part*/) { return true; },
            /*GridLayout const&*/ layout);
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
