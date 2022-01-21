#include "mkn/kul/log.hpp"

#include <thread>
#include <cstdlib>

#include "bench/core/bench.h"

#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"
#include "core/utilities/thread_pool.hpp"

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

    auto static constexpr range_maker(std::vector<IonPatchView>& views, std::uint32_t n_threads)
    {
        return PHARE::core::updater_ranges_per_thread(views, n_threads);
    }

public:
    void SetUp(::benchmark::State const& state) override
    {
        auto constexpr disperse = false;
        std::size_t n_threads   = state.threads();
        std::uint32_t n_parts   = state.range(0);
        std::uint32_t cells     = state.range(1);

        if (state.thread_index() == 0)
        {
            patches = core::generate(
                [&](std::size_t /*i*/) {
                    return PatchState::make_unique(
                        n_parts, PHARE::core::ConstArray<std::uint32_t, dim>(0),
                        PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1), disperse);
                },
                n_threads);

            views = core::generate(
                [&](auto& patch) -> IonPatchView {
                    return {patch->layout, patch->EM.view(), IonPopView::make_shared(patch->ions)};
                },
                patches);

            units = range_maker(views, state.threads());
        }
    }

    void TearDown(::benchmark::State const& state) override
    {
        if (state.thread_index() == 0) {}
    }

    void run(::benchmark::State& state)
    {
        Interpolator interpolator;

        for (auto _ : state)
        {
            for (auto& particles : units[state.thread_index()])
            {
                auto& pop    = *particles.domain.view;
                auto& layout = *particles.domain.layout;

                interpolator.template particleToMesh<version>(particles.domain, pop.density,
                                                              pop.flux, layout);
            }
        }
    }

    std::vector<std::unique_ptr<PatchState>> patches;
    std::vector<IonPatchView> views;
    decltype(range_maker(views, 0)) units;
};

BENCHMARK_TEMPLATE_DEFINE_F(Part2MeshFixture, _2_1_0, 2, 1, 0)(benchmark::State& state)
{
    run(state);
}
BENCHMARK_REGISTER_F(Part2MeshFixture, _2_1_0)
    ->Unit(benchmark::kNanosecond)
    ->Threads(1)
    ->Threads(2)
    // ->Threads(3)
    // ->Threads(4)
    ->Threads(5)
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
//     ->Unit(benchmark::kNanosecond)
//     ->Threads(1)
//     ->Threads(2)
//     ->Threads(3)
//     ->Threads(4)
//     ->Threads(5)
//     ->Threads(6)
//     ->Threads(7)
//     ->Threads(8)
//     ->Threads(9)
//     ->Threads(10)
//     ->Threads(15)
//     ->Threads(20)
//     ->Threads(25)
//     ->Threads(30)
//     ->UseRealTime()
//     ->ArgsProduct({TOT, CELLS});


// BENCHMARK_TEMPLATE_DEFINE_F(Part2MeshFixture, _2_1_2, 2, 1, 2)(benchmark::State& state)
// {
//     run(state);
// }
// BENCHMARK_REGISTER_F(Part2MeshFixture, _2_1_2)
//     ->Unit(benchmark::kNanosecond)
//     ->Threads(1)
//     ->Threads(2)
//     ->Threads(3)
//     ->Threads(4)
//     ->Threads(5)
//     ->Threads(6)
//     ->Threads(7)
//     ->Threads(8)
//     ->Threads(9)
//     ->Threads(10)
//     ->Threads(15)
//     ->Threads(20)
//     ->Threads(25)
//     ->Threads(30)
//     ->UseRealTime()
//     ->ArgsProduct({TOT, CELLS});



} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
