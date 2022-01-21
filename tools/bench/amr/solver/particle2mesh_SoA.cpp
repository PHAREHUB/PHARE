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

static std::vector<std::int64_t> TOT{500000, 1000000};
static std::vector<std::int64_t> CELLS{400};

template<std::size_t dim, std::size_t interp, bool sorted>
struct SolveFixture : public benchmark::Fixture
{
    auto constexpr atomic    = false;
    using Interpolator       = PHARE::core::Interpolator<dim, interp, atomic>;
    using PHARE_Types        = PHARE::core::PHARE_Types<dim, interp>;
    using ParticleArray_t    = typename PHARE_Types::ParticleSoA_t;
    using GridLayout_t       = typename PHARE_Types::GridLayout_t;
    using VecField_t         = typename PHARE_Types::VecField_t;
    using Electromag_t       = typename PHARE_Types::Electromag_t;
    using IonPopView         = core::IonPopulationView<ParticleArray_t, VecField_t, GridLayout_t>;
    using RangeSynchrotron_t = PHARE::core::RangeSynchrotron<ParticleArray_t>;
    using IonPatchView       = std::tuple<GridLayout_t, typename Electromag_t::view_t,
                                    std::vector<std::shared_ptr<IonPopView>>>;
    using PatchState =
        typename PHARE::core::bench::HybridPatch<GridLayout_t, ParticleArray_t>::State;


public:
    void SetUp(::benchmark::State const& state) override
    {
        std::uint32_t n_parts = state.range(0);
        std::uint32_t cells   = state.range(1);
        patches.emplace_back(
            PatchState::make_unique(n_parts, PHARE::core::ConstArray<std::uint32_t, dim>(0),
                                    PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1)));
        if constexpr (sorted)
            for (auto& patch : patches)
                for (auto& pop : patch->ions)
                    std::sort(pop.domainParticles());

        for (auto& patch : patches)
            ion_patch_views.emplace_back(patch->layout, patch->EM.view(),
                                         IonPopView::make_shared(patch->ions));
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

template<std::size_t dim, std::size_t interp, bool sorted>
void SolveFixture<dim, interp, sorted>::solve(::benchmark::State& state)
{
    assert(patches.size());

    Interpolator interpolator;
    auto units = PHARE::core::updater_ranges_per_thread(ion_patch_views);

    for (auto& particles : units[0])
    {
        auto& pop    = *particles.domain.view;
        auto& layout = *particles.domain.layout;
        auto& em     = *particles.domain.em;

        for (auto _ : state)
            interpolator.particleToMesh(particles.domain, pop.density, pop.flux, layout);
    }
}

BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _3_1_0_push, 3, 1, 0)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _3_1_0_push)
    ->Unit(benchmark::kNanosecond)
    ->ArgsProduct({TOT, CELLS});


BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _3_1_1_push, 3, 1, 1)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _3_1_1_push)
    ->Unit(benchmark::kNanosecond)
    ->ArgsProduct({TOT, CELLS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}