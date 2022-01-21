#include "mkn/kul/log.hpp"

#include <thread>
#include <cstdlib>

#include "bench/core/bench.h"

#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"
#include "core/utilities/thread_pool.hpp"

namespace PHARE::amr::bench
{
static std::vector<std::int64_t> THREADS{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12};
static std::vector<std::int64_t> TOT{500000, 1000000, 10000000};
static std::vector<std::int64_t> CELLS{200, 250, 300, 350, 400};

template<std::size_t dim, std::size_t interp, bool sorted>
struct SolveFixture : public benchmark::Fixture
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
    void SetUp(::benchmark::State const& state) override
    {
        auto constexpr disperse = false;
        std::uint32_t n_threads = state.range(0);
        std::uint32_t n_parts   = state.range(1);
        std::uint32_t cells     = state.range(2);
        patches.emplace_back(PatchState::make_unique(
            n_parts, PHARE::core::ConstArray<std::uint32_t, dim>(0),
            PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1), disperse));

        for (auto& patch : patches)
        {
            if constexpr (sorted)
                for (auto& pop : patch->ions)
                    std::sort(pop.domainParticles());

            ion_patch_views.emplace_back(patch->layout, patch->EM.view(),
                                         IonPopView::make_shared(patch->ions));
        }
    }

    void TearDown(::benchmark::State const& /*state*/) override
    {
        ion_patch_views.clear();
        patches.clear();
    }

    void solve(::benchmark::State&);

    std::vector<std::unique_ptr<PatchState>> patches;
    std::vector<IonPatchView> ion_patch_views;
};

template<std::size_t dim, std::size_t interp, bool sorted>
void SolveFixture<dim, interp, sorted>::solve(::benchmark::State& state)
{
    assert(patches.size());

    Interpolator interpolator;

    std::uint32_t n_threads = state.range(0);
    auto units              = PHARE::core::updater_ranges_per_thread(ion_patch_views, n_threads);

    auto thread_fn = [&](std::uint16_t idx) {
        for (auto& particles : units[idx])
        {
            auto& pop    = *particles.domain.view;
            auto& layout = *particles.domain.layout;
            auto& em     = *particles.domain.em;

            interpolator.template particleToMesh<1>(particles.domain, pop.density, pop.flux,
                                                    layout);
        }
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

BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_0_push, 2, 1, 0)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_0_push)
    ->Unit(benchmark::kNanosecond)
    ->ArgsProduct({THREADS, TOT, CELLS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}