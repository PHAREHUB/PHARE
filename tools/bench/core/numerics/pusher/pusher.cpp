#include "push_bench.hpp"

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    constexpr std::uint32_t cells   = 65;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types       = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t      = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Electromag_t      = PHARE::core::bench::Electromag<GridLayout_t>;
    using Ions_t            = typename PHARE_Types::Ions_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using Particle_t        = typename ParticleArray::value_type;
    using ParticleRange     = PHARE::core::IndexRange<ParticleArray>;

    using BorisPusher_t = PHARE::core::BorisPusher<dim, ParticleRange, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;


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

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    Interpolator interpolator;
    auto const no_op = [](auto& particleRange) { return particleRange; };
    while (state.KeepRunningBatch(1))
    {
        pusher.move(
            /*ParticleRange const&*/ rangeIn,  //
            /*ParticleRange&      */ rangeOut, //
            /*Electromag const&   */ em,       //
            /*double mass         */ 1,
            /*Interpolator&*/ interpolator,
            /*GridLayout const&*/ layout, //
            no_op, no_op);
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
