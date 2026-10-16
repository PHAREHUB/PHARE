
#include "core/numerics/interpolator/interpolator.hpp"

#include "simulator/simulator_def.hpp"

#include "tools/bench/core/bench.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "benchmark/benchmark.h"

template<std::size_t dim, std::size_t interp>
void interpolate(benchmark::State& state)
{
    constexpr std::uint32_t cells = 30;
    constexpr std::uint32_t ppc   = 100;
    auto static constexpr opts    = PHARE::SimOpts{dim, interp};

    using PHARE_Types   = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using ParticleArray = PHARE_Types::Hybrid::ParticleArray_t;
    using Grid_t        = PHARE_Types::Hybrid::Grid_t;
    using Hybrid_t      = PHARE_Types::Hybrid;

    auto constexpr static field_opts = PHARE::core::TensorFieldOptions<Hybrid_t>{};

    GridLayout_t layout{cells};
    PHARE::core::Interpolator<dim, interp> interpolator;
    ParticleArray particles{layout.AMRBox()};
    add_particles(particles, layout.AMRBox(), ppc);

    PHARE::core::UsableElectromag<field_opts> em{layout};
    PHARE::core::UsableVecField<field_opts> flux{"F", layout,
                                                 PHARE::core::HybridQuantity::Vector::V};
    Grid_t particleDensity{"particleDensity", PHARE::core::HybridQuantity::Scalar::rho,
                           layout.allocSize(PHARE::core::HybridQuantity::Scalar::rho)};
    Grid_t chargeDensity{"chargeDensity", PHARE::core::HybridQuantity::Scalar::rho,
                         layout.allocSize(PHARE::core::HybridQuantity::Scalar::rho)};

    delta_disperse(particles);

    while (state.KeepRunning())
    {
        // meshToParticle
        for (std::size_t i = 0; i < particles.size(); ++i)
            interpolator(particles, em, layout, i);

        // particleToMesh
        interpolator(particles, particleDensity, chargeDensity, flux, layout);
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
