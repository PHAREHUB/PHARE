#include "bench_interpolator.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"

namespace PHARE::core
{

template<std::size_t dim>
using AoSGPUParticleArray
    = ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::AoS, StorageMode::VECTOR,
                                                PHARE::AllocatorMode::GPU_UNIFIED>>;

template<std::size_t dim>
using SoAGPUParticleArray
    = ParticleArray<dim, ParticleArrayInternals<dim, LayoutMode::SoA, StorageMode::VECTOR,
                                                PHARE::AllocatorMode::GPU_UNIFIED>>;

} // namespace PHARE::core

template<std::size_t dim, std::size_t interp, bool sort>
void bench_cpu_soa(benchmark::State& state)
{
    interpolate<dim, interp, PHARE::core::SoAParticleArray<dim>, sort>(state);
}

template<std::size_t dim, std::size_t interp, bool sort>
void bench_cpu_aos(benchmark::State& state)
{
    interpolate<dim, interp, PHARE::core::AoSParticleArray<dim>, sort>(state);
}

template<std::size_t dim, std::size_t interp, bool sort>
void bench_gpu_soa(benchmark::State& state)
{
    interpolate<dim, interp, PHARE::core::SoAGPUParticleArray<dim>, sort>(state);
}

template<std::size_t dim, std::size_t interp, bool sort>
void bench_gpu_aos(benchmark::State& state)
{
    interpolate<dim, interp, PHARE::core::AoSGPUParticleArray<dim>, sort>(state);
}

BENCHMARK_TEMPLATE(bench_cpu_aos, 3, 3, false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(bench_gpu_aos, 3, 3, false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(bench_cpu_soa, 3, 3, false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(bench_gpu_soa, 3, 3, false)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(bench_cpu_aos, 3, 3, true)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(bench_gpu_aos, 3, 3, true)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(bench_cpu_soa, 3, 3, true)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(bench_gpu_soa, 3, 3, true)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
