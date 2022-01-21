


#include <omp.h>
#include "bench/core/bench.h"

#include "core/numerics/ion_updater/ion_range.h"
#include "core/numerics/ion_updater/ion_updater.h"
#include "core/utilities/thread_pool.hpp"

namespace PHARE::amr::bench
{
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
        constexpr std::uint32_t cells = 1e2;
        constexpr std::uint32_t parts = 1e3;
        std::uint32_t n_patches       = 10;
        for (std::size_t i = 0; i < n_patches; ++i)
            patches.emplace_back(PatchState::make_unique(
                PHARE::core::ConstArray<std::uint32_t, dim>(0),
                PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1), parts));
        assert(patches.size());
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
    std::uint16_t n_threads = omp_get_max_threads();

    for (auto& patch : patches)
        ion_patch_views.emplace_back(patch->layout, patch->EM.view(),
                                     IonPopView::make_shared(patch->ions));

    auto units = PHARE::core::updater_ranges_per_thread(ion_patch_views, n_threads);

    auto synchrotron = std::make_shared<RangeSynchrotron_t>(n_threads); // syncs on destruct

    auto thread_fn = [&](std::uint16_t idx) {
        IonUpdater_t{"modified_boris", idx, synchrotron}.updatePopulations(
            units[idx], 1e-5, PHARE::core::UpdaterMode::domain_only);
    };

    for (auto _ : state)
    {
#pragma omp parallel for
        for (auto& unit : units)
        {
            std::uint16_t t_idx = omp_get_thread_num();
            IonUpdater_t{"modified_boris", t_idx, synchrotron}.updatePopulations(
                unit, 1e-5, PHARE::core::UpdaterMode::domain_only);
        }
    }
}


BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_1_push, 2, 1, 1)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_1_push)->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_0_push, 2, 1, 0)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_0_push)->Unit(benchmark::kNanosecond);

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}