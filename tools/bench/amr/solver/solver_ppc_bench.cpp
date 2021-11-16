
#include "bench/amr/solver/solver_ppc_bench.h"

#include <atomic>
#include <thread>

#include "core/numerics/ohm/ohm.h"
#include "core/numerics/ampere/ampere.h"
#include "core/numerics/faraday/faraday.h"
#include "core/numerics/ion_updater/ion_updater.h"


#ifndef _PHARE_BENCH_THREAD_START_
#define _PHARE_BENCH_THREAD_START_ 1
#endif

#ifndef _PHARE_BENCH_THREADS_END_
#define _PHARE_BENCH_THREADS_END_ 20
#endif


namespace PHARE::amr::bench
{
std::size_t static constexpr thread_start = _PHARE_BENCH_THREAD_START_;
std::size_t static constexpr thread_end   = _PHARE_BENCH_THREADS_END_;

static auto threads()
{
    std::vector<std::int64_t> vec;
    for (std::size_t i = thread_start; i <= thread_end; ++i)
        vec.emplace_back(i);
    return vec;
}
static std::vector<std::int64_t> THREADS{threads()};

template<std::size_t dim, std::size_t interp, std::size_t op = 1>
struct SolveFixture : public benchmark::Fixture
{
    using PHARE_Types  = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t = typename PHARE_Types::GridLayout_t;
    using Ions_t       = typename PHARE_Types::Ions_t;
    using PatchState   = typename PHARE::core::bench::HybridPatch<GridLayout_t>::State;
    using Electromag_t = PHARE::core::bench::Electromag<GridLayout_t>;
    using IonUpdater   = PHARE::core::IonUpdater<Ions_t, Electromag_t, GridLayout_t>;

public:
    void SetUp(::benchmark::State const& state) override
    {
        constexpr std::uint32_t cells = 1e3;
        constexpr std::uint32_t parts = 1e7;
        assert(state.range(0) > 0);
        std::uint32_t n_threads = static_cast<std::uint32_t>(state.range(0));
        assert(n_threads);
        for (std::size_t i = 0; i < n_threads; ++i)
            patches.emplace_back(PatchState::make_unique(
                PHARE::core::ConstArray<std::uint32_t, dim>(0),
                PHARE::core::ConstArray<std::uint32_t, dim>(cells - 1), parts));
        assert(patches.size());
    }

    void TearDown(::benchmark::State const& /*state*/) override { patches.clear(); }

    void solve(::benchmark::State&);

    std::vector<std::unique_ptr<PatchState>> patches;
};

template<std::size_t dim, std::size_t interp, std::size_t op>
void SolveFixture<dim, interp, op>::solve(::benchmark::State& state)
{
    using Faraday  = PHARE::core::Faraday<GridLayout_t>;
    using Ampere   = PHARE::core::Ampere<GridLayout_t>;
    using Ohm      = PHARE::core::Ohm<GridLayout_t>;
    using VecField = PHARE::core::bench::_VF_<GridLayout_t>;

    assert(patches.size());
    auto n_threads = state.range(0);
    assert(n_threads);
    auto units = PHARE::amr::bench::decompose(patches, n_threads);
    assert(units.size());
    assert(units.size() == n_threads);

    // corrector
    auto per_thread = [&](std::size_t idx) {
        //        KLOG(INF) << idx;
        IonUpdater ionUpdater_{"modified_boris"};
        for (auto& view : units[idx])
        {
            [[maybe_unused]] VecField& E  = view.EM->E;
            [[maybe_unused]] VecField& B  = view.EM->B;
            [[maybe_unused]] VecField& J  = (*view.J);
            [[maybe_unused]] auto& rho    = (*view.rho);
            [[maybe_unused]] auto& layout = view.layout;

            if constexpr (op == 1)
                for (std ::size_t i = 0; i < 1e2; ++i)
                {
                    Faraday::op(layout, B, E, B, .01, view.B_boxes());
                    Ampere::op(layout, B, J, view.J_boxes());
                    Ohm::op(layout, rho, E, rho, B, J, E, view.E_boxes());
                }

            if constexpr (op == 2)
                for (std ::size_t i = 0; i < 1; ++i)
                    ionUpdater_.updatePopulations(*view.ions, *view.EM, layout, 1e-5,
                                                  PHARE::core::UpdaterMode::domain_only);
        }
    };


    for (auto _ : state)
    {
        auto threads = PHARE::core::generate(
            [&](auto i) { return std::thread{[&, i]() { per_thread(i); }}; }, n_threads - 1);
        per_thread(n_threads - 1);

        for (auto& thread : threads)
            if (thread.joinable())
                thread.join();
    }
}

// BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1, 2, 1)(benchmark::State& state)
// {
//     solve(state);
// }
// BENCHMARK_REGISTER_F(SolveFixture, _2_1)->Unit(benchmark::kNanosecond)->ArgsProduct({THREADS});
// other interps are not really useful


BENCHMARK_TEMPLATE_DEFINE_F(SolveFixture, _2_1_push, 2, 1, 2)(benchmark::State& state)
{
    solve(state);
}
BENCHMARK_REGISTER_F(SolveFixture, _2_1_push)->Unit(benchmark::kNanosecond)->ArgsProduct({THREADS});
} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
