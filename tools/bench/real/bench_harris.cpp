
#include "bench/real/bench_harris.h"

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

template<std::size_t dim, std::size_t interp, std::size_t nbRefinePart>
struct SimulatorFixture : public benchmark::Fixture
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
        if (!hierarchy)
        {
            auto& dict = PHARE::real::bench::harris::createDict();
            hierarchy  = PHARE::amr::Hierarchy::make();
            simulator
                = std::make_unique<typename decltype(simulator)::element_type>(dict, hierarchy);
            simulator->initialize();
            PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        }
    }

    void TearDown(::benchmark::State const& /*state*/) override
    {
        // hierarchy.reset();
        // simulator.reset();
        // PHARE::SamraiLifeCycle::reset();
    }

    void run(::benchmark::State&);

    std::shared_ptr<Hierarchy> hierarchy;
    std::unique_ptr<PHARE::Simulator<dim, interp, nbRefinePart>> simulator;
};

template<std::size_t dim, std::size_t interp, std::size_t op>
void SimulatorFixture<dim, interp, op>::run(::benchmark::State& state)
{
    KLOG(INF);
    for (auto _ : state)
    {
        // while (simulator->currentTime() < simulator->endTime())
        // {
        //     // simulator->dump(simulator->currentTime(), simulator->timeStep());
        simulator->advance(simulator->timeStep());
        // }
    }
}

BENCHMARK_TEMPLATE_DEFINE_F(SimulatorFixture, _2_1_advance, 2, 1, 4)(benchmark::State& state)
{
    run(state);
}
BENCHMARK_REGISTER_F(SimulatorFixture, _2_1_advance)->Unit(benchmark::kNanosecond);
// ->ArgsProduct({THREADS});

} // namespace PHARE::amr::bench

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
