#include "mkn/kul/log.hpp"

#include <thread>
#include <cstdlib>

#include "bench/core/bench.h"

#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"

namespace PHARE::amr::bench
{
static std::vector<std::int64_t> TOT{/*100000,*/ 5000000};
static std::vector<std::int64_t> CELLS{/*100, */ 400};

template<std::size_t dim, std::size_t interp, std::size_t version>
struct Part2MeshFixture : public benchmark::Fixture
{
    auto constexpr atomic = false;
    using Interpolator    = PHARE::core::Interpolator<dim, interp, atomic>;
    using PHARE_Types     = PHARE::core::PHARE_Types<dim, interp>;
    using ParticleArray_t = typename PHARE_Types::ParticleSoA_t;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using VecField_t      = typename PHARE_Types::VecField_t;
    using Electromag_t    = typename PHARE_Types::Electromag_t;
    using IonPopView      = core::IonPopulationView<ParticleArray_t, VecField_t, GridLayout_t>;
    using IonPatchView    = std::tuple<GridLayout_t, typename Electromag_t::view_t,
                                    std::vector<std::shared_ptr<IonPopView>>>;
    using PatchState =
        typename PHARE::core::bench::HybridPatch<GridLayout_t, ParticleArray_t>::State;


public:
    void SetUp(::benchmark::State const& state) override {}

    void TearDown(::benchmark::State const& state) override {}

    void run(::benchmark::State& state)
    {
        auto constexpr disperse = true;
        std::size_t n_threads   = 1; // state.threads;
        std::uint32_t n_parts   = state.range(0);
        std::uint32_t cells     = state.range(1);

        auto patches = core::generate(
            [&](std::size_t /*i*/) {
                return PatchState::make_unique(
                    n_parts, PHARE::core::ConstArray<std::uint32_t, dim>(0),
                    PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1), disperse);
            },
            n_threads);

        for (auto& patch : patches)
            for (auto& pop : patch->ions)
                std::sort(pop.domainParticles());

        auto views = core::generate(
            [&](auto& patch) -> IonPatchView {
                return {patch->layout, patch->EM.view(), IonPopView::make_shared(patch->ions)};
            },
            patches);

        Interpolator interpolator;
        auto units = PHARE::core::updater_ranges_per_thread(views);
        assert(units.size() == 1); // one per thread

        for (auto _ : state)
        {
            for (auto& particles : units[0])
            {
                auto& pop    = *particles.domain.view;
                auto& layout = *particles.domain.layout;

                interpolator.template particleToMesh<version>(particles.domain, pop.density,
                                                              pop.flux, layout);
            }
        }
    }
};

BENCHMARK_TEMPLATE_DEFINE_F(Part2MeshFixture, _2_1_0, 2, 1, 0)(benchmark::State& state)
{
    run(state);
}
BENCHMARK_REGISTER_F(Part2MeshFixture, _2_1_0)
    ->Unit(benchmark::kNanosecond)
    // ->Threads(1)
    // ->Threads(2)
    // ->Threads(3)
    // ->Threads(4)
    // ->Threads(5)
    // ->Threads(6)
    // ->Threads(7)
    // ->Threads(8)
    // ->Threads(9)
    // ->Threads(10)
    // ->Threads(15)
    ->Threads(20)
    // ->Threads(25)
    // ->Threads(30)
    ->ArgsProduct({TOT, CELLS});


// BENCHMARK_TEMPLATE_DEFINE_F(Part2MeshFixture, _2_1_1, 2, 1, 1)(benchmark::State& state)
// {
//     run(state);
// }
// BENCHMARK_REGISTER_F(Part2MeshFixture, _2_1_1)
//     ->MeasureProcessCPUTime()
//     ->Unit(benchmark::kNanosecond)
//     ->Threads(1)
//     ->Threads(2)
//     // ->Threads(3)
//     // ->Threads(4)
//     ->Threads(5)
//     // ->Threads(6)
//     // ->Threads(7)
//     // ->Threads(8)
//     // ->Threads(9)
//     // ->Threads(10)
//     // ->Threads(15)
//     ->Threads(20)
//     // ->Threads(25)
//     // ->Threads(30)
//     ->ArgsProduct({TOT, CELLS});


// BENCHMARK_TEMPLATE_DEFINE_F(Part2MeshFixture, _2_1_2, 2, 1, 2)(benchmark::State& state)
// {
//     run(state);
// }
// BENCHMARK_REGISTER_F(Part2MeshFixture, _2_1_2)
//     ->MeasureProcessCPUTime()
//     ->Unit(benchmark::kNanosecond)
//     ->Threads(1)
//     ->Threads(2)
//     // ->Threads(3)
//     // ->Threads(4)
//     ->Threads(5)
//     // ->Threads(6)
//     // ->Threads(7)
//     // ->Threads(8)
//     // ->Threads(9)
//     // ->Threads(10)
//     // ->Threads(15)
//     ->Threads(20)
//     // ->Threads(25)
//     // ->Threads(30)
//     ->ArgsProduct({TOT, CELLS});


} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
