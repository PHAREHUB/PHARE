
#include "tools/bench/core/bench.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

template<std::size_t dim, std::size_t interp>
void interpolate(benchmark::State& state)
{
    constexpr std::uint32_t cells   = 30;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types   = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t  = typename PHARE_Types::GridLayout_t;
    using ParticleArray = typename PHARE_Types::ParticleArray_t;

    PHARE::core::Interpolator<dim, interp> interpolator;
    ParticleArray particles{n_parts, PHARE::core::bench::particle<dim>()};
    GridLayout_t layout{PHARE::core::ConstArray<double, dim>(1.0 / cells),
                        PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                        PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}};
    PHARE::core::bench::Electromag<GridLayout_t> em{layout};
    PHARE::core::bench::Flux<GridLayout_t> flux{layout};
    auto rho = PHARE::core::bench::rho(layout);

    PHARE::core::bench::disperse(particles, 0, cells - 1);

    while (state.KeepRunning())
    {
        // meshToParticle
        interpolator(particles, em, layout);

        // particleToMesh
        interpolator(particles, rho, flux, layout);
    }
}

BENCHMARK_TEMPLATE(interpolate, 1, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(interpolate, 1, 2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(interpolate, 1, 3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(interpolate, 2, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(interpolate, 2, 2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(interpolate, 2, 3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(interpolate, 3, 1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(interpolate, 3, 2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(interpolate, 3, 3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
