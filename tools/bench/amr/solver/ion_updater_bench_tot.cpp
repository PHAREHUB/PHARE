

#include "mkn/kul/log.hpp"

#include "bench/core/bench.h"

#include <thread>

#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"
#include "core/utilities/thread_pool.hpp"



namespace PHARE::amr::bench
{
static std::vector<std::int64_t> THREADS{1, 2, 3, 4, 5, 10, 12, 15, 20};
static std::vector<std::int64_t> TOT{10000, 20000, 30000, 50000, 100000};
static std::vector<std::int64_t> CELLS{1, 2, 3, 4, 5, 10, 12, 15, 20, 30, 50};

template<std::size_t dim, std::size_t interp, bool atomic = true>
struct SolveFixture : public benchmark::Fixture
{
    using PHARE_Types = PHARE::core::PHARE_Types<dim, interp>;

    using ParticleArray_t = typename PHARE_Types::ParticleArray_t;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using VecField_t      = typename PHARE_Types::VecField_t;
    using PatchState      = typename PHARE::core::bench::HybridPatch<GridLayout_t>::State;
    using Electromag_t    = typename PHARE_Types::Electromag_t;

    using IonPopView = core::IonPopulationView<ParticleArray_t, VecField_t, GridLayout_t>;
    using IonUpdater_t
        = core::IonUpdater<typename Electromag_t::view_t, ParticleArray_t, GridLayout_t, atomic>;
    using RangeSynchrotron_t = PHARE::core::RangeSynchrotron<ParticleArray_t>;
    using IonPatchView       = std::tuple<GridLayout_t, typename Electromag_t::view_t,
                                    std::vector<std::shared_ptr<IonPopView>>>;


public:
    void SetUp(::benchmark::State const& state) override
    {
        std::uint32_t n_threads = state.range(0);
        std::uint32_t n_parts   = state.range(1);
        std::uint32_t cells     = state.range(2);
        patches.emplace_back(
            PatchState::make_unique(n_parts, PHARE::core::ConstArray<std::uint32_t, dim>(0),
                                    PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1)));
    }

    void TearDown(::benchmark::State const& /*state*/) override
    {
        patches.clear();
        ion_patch_views.clear();
    }

    void solve(::benchmark::State&);

    std::vector<std::unique_ptr<PatchState>> patches;
    std::vector<IonPatchView> ion_patch_views;
};

template<std::size_t dim, std::size_t interp, bool atomic>
void SolveFixture<dim, interp, atomic>::solve(::benchmark::State& state)
{
    assert(patches.size());
    auto n_threads_ = state.range(0);
    core::abort_if(n_threads_ < 1);
    std::uint16_t n_threads = n_threads_;

    for (auto& patch : patches)
        ion_patch_views.emplace_back(patch->layout, patch->EM.view(),
                                     IonPopView::make_shared(patch->ions));

    auto units = PHARE::core::updater_ranges_per_thread(ion_patch_views, n_threads);

    auto synchrotron = std::make_shared<RangeSynchrotron_t>(n_threads); // syncs on destruct

    auto thread_fn = [&](std::uint16_t idx) {
        IonUpdater_t{"modified_boris", idx, synchrotron}.updatePopulations(
            units[idx], 1e-5, PHARE::core::UpdaterMode::domain_only);
    };

    thread_pool thread_pool_(n_threads - 1);

    for (auto _ : state)
    {
        for (std::size_t i = 1; i < n_threads; ++i)
            thread_pool_.submit([&, i]() { thread_fn(i); });

        thread_fn(0);
        thread_pool_.wait_for_tasks();
    }
}


BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_1_push, 2, 1, 1)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_1_push)
    ->Unit(benchmark::kNanosecond)
    ->ArgsProduct({THREADS, TOT, CELLS});

BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_0_push, 2, 1, 0)(benchmark::State& state)
{
    solve(state);
}
// BENCHMARK_REGISTER_F(SolveFixture, _2_1_0_push)
//     ->Unit(benchmark::kNanosecond)
//     ->ArgsProduct({THREADS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}