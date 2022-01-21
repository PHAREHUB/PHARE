// #include "mkn/kul/log.hpp"

#include <thread>
#include <cstdlib>

#include "bench/core/bench.h"

#include "core/data/particles/particle_array_soa.h"

namespace PHARE::amr::bench
{
using namespace PHARE::core;

static std::vector<std::int64_t> TOT{500000, 1000000};
static std::vector<std::int64_t> CELLS{100, 200, 300, 400, 500};

template<std::size_t dim, std::size_t interp>
struct SortFixture : public benchmark::Fixture
{
    using PHARE_Types = core::PHARE_Types<dim, interp>;

    using ParticleArray_t = core::ParticleArray_SOA<dim>;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using PatchState      = typename core::bench::HybridPatch<GridLayout_t, ParticleArray_t>::State;

public:
    void SetUp(::benchmark::State const& state) override
    {
        std::uint32_t n_parts = state.range(0);
        std::uint32_t cells   = state.range(1);

        patch = PatchState::make_unique(n_parts, ConstArray<std::uint32_t, dim>(0),
                                        ConstArray<std::uint32_t, dim>(cells - 1));
    }

    void operator()(::benchmark::State& state)
    {
        for (auto& pop : patch->ions)
            for (auto _ : state)
                std::sort(pop.domainParticles());

        for (auto& pop : patch->ions)
        {
            auto prev = &pop.domainParticles()[0];
            for (std::size_t i = 1; i < pop.domainParticles().size(); ++i)
            {
                auto& next = pop.domainParticles()[i];

                abort_if(prev->iCell > next.iCell);

                prev = &next;
            }
        }
    }

    std::unique_ptr<PatchState> patch;
};

BENCHMARK_TEMPLATE_DEFINE_F(SortFixture, _3_1, 3, 1)(benchmark::State& state)
{
    (*this)(state);
}
BENCHMARK_REGISTER_F(SortFixture, _3_1)->Unit(benchmark::kNanosecond)->ArgsProduct({TOT, CELLS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}