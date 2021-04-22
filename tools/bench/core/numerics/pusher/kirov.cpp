
#include "benchmark/benchmark.h"
#include "pusher_bench.h"

using namespace PHARE::core::bench;

static std::size_t nThreads = 1;

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    constexpr std::uint32_t cells = 65;
    constexpr std::uint32_t parts = 1e8;

    using PHARE_Types       = PHARE::core::PHARE_Types<dim, interp>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Ions_t            = typename PHARE_Types::Ions_t;
    using Electromag_t      = typename PHARE_Types::Electromag_t;
    using GridLayout_t      = typename PHARE_Types::GridLayout_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using PartIterator      = typename ParticleArray::iterator;
    using ThreadPool        = ::EXT::ThreadPool;

    using KirovPusher_t = PHARE::core::KirovPusher<dim, PartIterator, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t, ThreadPool>;

    Interpolator interpolator;
    ParticleArray domainParticles{parts, particle<dim>(/*icell =*/34)};

    auto range    = PHARE::core::makeRange(domainParticles);
    auto meshSize = PHARE::core::ConstArray<double, dim>(1.0 / cells);
    auto nCells   = PHARE::core::ConstArray<std::uint32_t, dim>(cells);
    auto origin   = PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)};

    GridLayout_t layout{meshSize, nCells, origin};

    PHARE::core::bench::Electromag<GridLayout_t, VecField<dim>> electromag{layout};

    KirovPusher_t pusher{nThreads};
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    while (state.KeepRunning())
    {
        pusher.move(
            /*ParticleRange&*/ range,
            /*Electromag const&*/ electromag,
            /*double mass*/ 1,
            /*Interpolator&*/ interpolator,
            /*ParticleSelector const&*/ [](auto const& /*part*/) { return true; },
            /*GridLayout const&*/ layout);
    }
}

BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);



void set_env(char** envp)
{
    const std::string find = "PHARE_THREADS=";
    for (char** environ = envp; *environ != 0; environ++)
        if (std::string env(*environ); env.find(find) != std::string::npos)
        {
            std::stringstream sstream(env.substr(find.size()));
            sstream >> nThreads;
        }
    std::cout << __FILE__ << " " << __LINE__ << " WITH PHARE_THREADS=" << nThreads << std::endl;
}

int main(int argc, char** argv, char** envp)
{
    set_env(envp);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
