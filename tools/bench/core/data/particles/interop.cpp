
#include "benchmark/benchmark.h"

#include "tools/bench/core/bench.hpp"
#include "core/data/particles/particle.hpp"

constexpr std::size_t dim  = 3;
constexpr std::size_t size = 100000;

auto makeContiguous()
{
    PHARE::core::ContiguousParticles<dim> contiguous{size};
    for (std::size_t i = 0; i < size; i++)
    {
        auto view   = contiguous[i];
        view.weight = 1 + i;
        view.charge = 1 + i;
        view.iCell  = PHARE::core::ConstArray<int, dim>(i);
        view.delta  = PHARE::core::ConstArray<double, dim>(i + 1);
        view.v      = PHARE::core::ConstArray<double, 3>(view.weight + 2);
    }
    return contiguous;
}

void contiguousRange(benchmark::State& state)
{
    auto contiguous = makeContiguous();

    while (state.KeepRunning())
        for (auto const& view : contiguous)
            benchmark::DoNotOptimize(view == std::copy(view));
}
BENCHMARK(contiguousRange)->Unit(benchmark::kMicrosecond);


void contiguousBracketOperator(benchmark::State& state)
{
    auto contiguous = makeContiguous();

    while (state.KeepRunning())
        for (std::size_t i = 0; i < contiguous.size(); i++)
        {
            auto view = contiguous[i];
            benchmark::DoNotOptimize(view == std::copy(view));
        }
}
BENCHMARK(contiguousBracketOperator)->Unit(benchmark::kMicrosecond);


int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
