
#include "bench/core/bench.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"



template<typename ParticleArray, std::size_t impl_ = 0>
struct Fixture_t : public benchmark::Fixture
{
    auto constexpr static dim            = ParticleArray::dimension;
    auto constexpr static impl           = impl_;
    std::uint32_t constexpr static cells = 50;
    std::uint32_t constexpr static ppc   = 100;

    using PHARE_Types     = PHARE::core::PHARE_Types<dim, 1>;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using ParticleArray_t = ParticleArray;

    void SetUp(::benchmark::State& state)
    {
        particles = std::make_unique<ParticleArray>();
        add_particles_in(*particles, layout.AMRBox(), ppc);
        PHARE::core::disperse(*particles, 1333337);
    }

    void TearDown(::benchmark::State& state) { particles.reset(); }

    std::unique_ptr<ParticleArray_t> particles;
    TestGridLayout<GridLayout_t> layout{cells};
};


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

BENCHMARK_TEMPLATE_DEFINE_F(Fixture_t, bench_aos_cpu, PHARE::core::AoSParticleArray<3>)
(benchmark::State& st)
{
    while (st.KeepRunningBatch(1))
        PHARE::core::ParticleArraySorter<ParticleArray_t, impl>{*particles, layout.AMRBox()}();
}
BENCHMARK_REGISTER_F(Fixture_t, bench_aos_cpu)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE_DEFINE_F(Fixture_t, bench_soa_cpu, PHARE::core::SoAParticleArray<3>)
(benchmark::State& st)
{
    while (st.KeepRunningBatch(1))
        PHARE::core::ParticleArraySorter<ParticleArray_t, impl>{*particles, layout.AMRBox()}();
}
BENCHMARK_REGISTER_F(Fixture_t, bench_soa_cpu)->Unit(benchmark::kMicrosecond);


BENCHMARK_TEMPLATE_DEFINE_F(Fixture_t, bench_aos_gpu, PHARE::core::AoSGPUParticleArray<3>)
(benchmark::State& st)
{
    while (st.KeepRunningBatch(1))
        PHARE::core::ParticleArraySorter<ParticleArray_t, impl>{*particles, layout.AMRBox()}();
}
BENCHMARK_REGISTER_F(Fixture_t, bench_aos_gpu)->Unit(benchmark::kMicrosecond);


BENCHMARK_TEMPLATE_DEFINE_F(Fixture_t, bench_soa_gpu, PHARE::core::SoAGPUParticleArray<3>)
(benchmark::State& st)
{
    while (st.KeepRunningBatch(1))
        PHARE::core::ParticleArraySorter<ParticleArray_t, impl>{*particles, layout.AMRBox()}();
}
BENCHMARK_REGISTER_F(Fixture_t, bench_soa_gpu)->Unit(benchmark::kMicrosecond);



int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
