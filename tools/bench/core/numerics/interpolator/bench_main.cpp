
#include "tools/bench/core/bench.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

template<std::size_t dim, std::size_t interp>
void interpolate(benchmark::State& state)
{
    constexpr std::uint32_t cells   = 30;
    constexpr std::uint32_t n_parts = 1e7;
    auto static constexpr opts      = PHARE::SimOpts<>{dim, interp};

    using PHARE_Types   = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using ParticleArray = PHARE_Types::ParticleArray_t;
    using Grid_t        = PHARE_Types::Grid_t;

    GridLayout_t layout{cells};
    PHARE::core::Interpolator<dim, interp> interpolator;
    ParticleArray particles{layout.AMRBox()};
    particles.vector().resize(n_parts, PHARE::core::bench::particle<dim>());
    PHARE::core::UsableElectromag<dim> em{layout};
    PHARE::core::UsableVecField<dim> flux{"F", layout, PHARE::core::HybridQuantity::Vector::V};
    Grid_t particleDensity{"particleDensity", PHARE::core::HybridQuantity::Scalar::rho,
                           layout.allocSize(PHARE::core::HybridQuantity::Scalar::rho)};
    Grid_t chargeDensity{"chargeDensity", PHARE::core::HybridQuantity::Scalar::rho,
                         layout.allocSize(PHARE::core::HybridQuantity::Scalar::rho)};

    PHARE::core::bench::disperse(particles, 0, cells - 1);

    while (state.KeepRunning())
    {
        // meshToParticle
        for (auto& particle : particles)
            interpolator(particle, em, layout);

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
