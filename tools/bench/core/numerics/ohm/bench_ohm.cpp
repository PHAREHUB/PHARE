#include "simulator/simulator_def.hpp"

#include "tools/bench/core/bench.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include <cmath>
#include <cstdint>

enum class OhmTerms { ideal, hyper_constant, hyper_spatial };

PHARE::core::OhmInfo ohm_info(OhmTerms const terms)
{
    using PHARE::core::HyperMode;
    if (terms == OhmTerms::ideal)
        return {/*eta=*/0., /*nu=*/0., HyperMode::constant};
    if (terms == OhmTerms::hyper_constant)
        return {/*eta=*/1., /*nu=*/0.01, HyperMode::constant};
    return {/*eta=*/1., /*nu=*/0.01, HyperMode::spatial};
}

// smooth, non-zero values everywhere (ghosts included): n divides in the pressure term
void fill(auto&& field, double const offset, double const phase)
{
    std::size_t i = 0;
    for (auto& v : field)
        v = offset + 0.5 * std::sin(0.01 * static_cast<double>(i++) + phase);
}

template<std::size_t dim, OhmTerms terms>
void ohm(benchmark::State& state)
{
    auto const cells           = static_cast<std::uint32_t>(state.range(0));
    auto static constexpr opts = PHARE::SimOpts{dim, 1};

    using PHARE_Types  = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using Grid_t       = PHARE_Types::Hybrid::Grid_t;
    using PHARE::core::HybridQuantity;

    GridLayout_t layout{cells};
    Grid_t n{"n", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)};
    Grid_t Pe{"Pe", HybridQuantity::Scalar::P, layout.allocSize(HybridQuantity::Scalar::P)};
    PHARE::core::UsableVecField<dim> Ve{"Ve", layout, HybridQuantity::Vector::V};
    PHARE::core::UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    PHARE::core::UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    PHARE::core::UsableVecField<dim> Enew{"Enew", layout, HybridQuantity::Vector::E};

    fill(n, 2., 0.);
    fill(Pe, 1., 0.5);
    for (std::size_t c = 0; c < 3; ++c)
    {
        fill(Ve[c], 0., 1. + c);
        fill(B[c], 1., 2. + c);
        fill(J[c], 0., 3. + c);
    }

    PHARE::core::Ohm<typename PHARE_Types::Hybrid::GridLayout_t> ohm{ohm_info(terms), layout};

    for (auto _ : state)
    {
        ohm(*&n, Ve, *&Pe, B, J, Enew);
        benchmark::DoNotOptimize(Enew[0].data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * 3 * std::pow(cells, dim));
}

#define PHARE_OHM_CASES(dim, small, large)                                                         \
    BENCHMARK_TEMPLATE(ohm, dim, OhmTerms::ideal)                                                  \
        ->Arg(small)                                                                               \
        ->Arg(large)                                                                               \
        ->Unit(benchmark::kMicrosecond);                                                           \
    BENCHMARK_TEMPLATE(ohm, dim, OhmTerms::hyper_constant)                                         \
        ->Arg(small)                                                                               \
        ->Arg(large)                                                                               \
        ->Unit(benchmark::kMicrosecond);                                                           \
    BENCHMARK_TEMPLATE(ohm, dim, OhmTerms::hyper_spatial)                                          \
        ->Arg(small)                                                                               \
        ->Arg(large)                                                                               \
        ->Unit(benchmark::kMicrosecond)

PHARE_OHM_CASES(2, 256, 1024);
PHARE_OHM_CASES(3, 32, 128);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
