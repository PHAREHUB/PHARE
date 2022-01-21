// #include "mkn/kul/log.hpp"

#include <thread>
#include <cstdlib>

#include "bench/core/bench.h"

#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"
#include "core/utilities/thread_pool.hpp"

namespace PHARE::amr::bench
{
auto cells()
{
    std::vector<std::int64_t> cells;
    for (std::size_t i = 4; i < 400; i += 2)
        cells.emplace_back(i);
    return cells;
}

static std::vector<std::int64_t> THREADS{1 /*, 9, 10, 11, 12, 13, 14, 15, 16, 20*/};
static std::vector<std::int64_t> TOT{/*100000,*/ 5000000};
static std::vector<std::int64_t> CELLS{/*100, */ 400};

template<std::size_t dim, std::size_t interp, bool atomic>
struct SolveFixture : public benchmark::Fixture
{
    using Interpolator    = PHARE::core::Interpolator<dim, interp, atomic>;
    using PHARE_Types     = PHARE::core::PHARE_Types<dim, interp>;
    using ParticleArray_t = typename PHARE_Types::ParticleArray_t;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using VecField_t      = typename PHARE_Types::VecField_t;
    using Electromag_t    = typename PHARE_Types::Electromag_t;
    using IonPopView      = core::IonPopulationView<ParticleArray_t, VecField_t, GridLayout_t>;
    using IonPatchView    = std::tuple<GridLayout_t, typename Electromag_t::view_t,
                                    std::vector<std::shared_ptr<IonPopView>>>;
    using PatchState      = typename PHARE::core::bench::HybridPatch<GridLayout_t>::State;


public:
    void SetUp(::benchmark::State const& state) override {}

    void TearDown(::benchmark::State const& /*state*/) override {}

    void solve(::benchmark::State&);
};

template<std::size_t dim, std::size_t interp, bool atomic>
void SolveFixture<dim, interp, atomic>::solve(::benchmark::State& state)
{
    std::vector<std::unique_ptr<PatchState>> patches;
    std::vector<IonPatchView> ion_patch_views;

    std::uint32_t n_threads = 1; // state.range(0);
    std::uint32_t n_parts   = state.range(1);
    std::uint32_t cells     = state.range(2);
    patches.emplace_back(
        PatchState::make_unique(n_parts, PHARE::core::ConstArray<std::uint32_t, dim>(0),
                                PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1)));

    assert(patches.size());
    // for (auto& patch : patches)
    //     for (auto& pop : patch->ions)
    //         std::sort(pop.domainParticles().begin(), pop.domainParticles().end(),
    //                   [](auto& x, auto& y) { return x.iCell > y.iCell; });

    for (auto& patch : patches)
        ion_patch_views.emplace_back(patch->layout, patch->EM.view(),
                                     IonPopView::make_shared(patch->ions));

    Interpolator interpolator;

    // std::uint32_t n_threads = state.range(0);
    auto units = PHARE::core::updater_ranges_per_thread(ion_patch_views, n_threads);

    auto thread_fn = [&](std::uint16_t idx) {
        for (auto& particles : units[idx])
        {
            auto& pop    = *particles.domain.view;
            auto& layout = *particles.domain.layout;
            auto& em     = *particles.domain.em;

            interpolator.particleToMesh(particles.domain, pop.density, pop.flux, layout);
        }
    };

    // thread_pool thread_pool_(n_threads - 1);

    for (auto _ : state)
    {
        // for (std::size_t i = 1; i < n_threads; ++i)
        //     thread_pool_.submit([&, i]() { thread_fn(i); });

        thread_fn(0);
        // thread_pool_.wait_for_tasks();
    }
}

BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_0_push, 2, 1, 0)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_0_push)
    ->Unit(benchmark::kNanosecond)
    ->Threads(1)
    ->Threads(30)
    ->UseRealTime()
    ->ArgsProduct({THREADS, TOT, CELLS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}