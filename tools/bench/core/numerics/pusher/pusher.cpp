
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"
#include "tools/bench/core/bench.hpp"

#include "benchmark/benchmark.h"

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    auto static constexpr opts = PHARE::SimOpts{dim, interp};
    using namespace PHARE::core;
    constexpr double mass         = 1;
    constexpr std::uint32_t cells = 10;
    constexpr std::uint32_t ppc   = 100;
    auto constexpr no_op          = [](auto& particleRange) { return particleRange; };

    using PHARE_Types       = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t      = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Electromag_t      = PHARE::core::UsableElectromag<GridLayout_t>;
    using Ions_t            = PHARE_Types::Ions_t;
    using ParticleArray_t   = Ions_t::particle_array_type;
    using ParticleRange     = PHARE::core::IndexRange<ParticleArray_t>;

    using BorisPusher_t = PHARE::core::BorisPusher<dim, ParticleRange, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;


    GridLayout_t layout{cells};
    Electromag_t em{layout};
    ParticleArray_t domainParticles{layout.AMRBox()};
    add_particles_in(domainParticles, layout.AMRBox(), ppc);
    delta_disperse(domainParticles, 13337);

    ParticleArray_t tmpDomain = domainParticles;
    auto rangeIn              = makeIndexRange(domainParticles);
    auto rangeOut             = makeIndexRange(tmpDomain);
    BorisPusher_t pusher{layout.meshSize(), .001};
    Interpolator interpolator;

    // while (state.KeepRunningBatch(1))
    pusher.move(rangeIn, rangeOut, em, mass, interpolator, layout, no_op, no_op);
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
