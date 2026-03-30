#define PHARE_UPDATER_IMPL 1

#include "push_bench.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"


#if PHARE_UPDATER_IMPL == 0

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    auto static constexpr opts      = PHARE::SimOpts{dim, interp};
    constexpr std::uint32_t cells   = 65;
    constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types       = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t      = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Electromag_t      = PHARE::core::UsableElectromag<dim>;
    using Ions_t            = PHARE_Types::Ions_t;
    using ParticleArray     = Ions_t::particle_array_type;
    using Particle_t        = ParticleArray::value_type;
    using ParticleRange     = PHARE::core::IndexRange<ParticleArray>;

    using BorisPusher_t = PHARE::core::BorisPusher<dim, ParticleRange, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;


    GridLayout_t layout{cells};
    Electromag_t em{layout};

    ParticleArray domainParticles{layout.AMRBox()};
    domainParticles.vector()
        = std::vector<Particle_t>(n_parts, PHARE::core::bench::particle<dim>());
    PHARE::core::bench::disperse(domainParticles, 0, cells - 1, 13337);
    // std::sort(domainParticles);

    ParticleArray tmpDomain{layout.AMRBox()};
    tmpDomain.vector() = domainParticles.vector();

    auto rangeIn  = PHARE::core::makeIndexRange(domainParticles);
    auto rangeOut = PHARE::core::makeIndexRange(tmpDomain);

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    Interpolator interpolator;
    auto const no_op = [](auto& particleRange) { return particleRange; };
    while (state.KeepRunningBatch(1))
    {
        pusher.move(
            /*ParticleRange const&*/ rangeIn,  //
            /*ParticleRange&      */ rangeOut, //
            /*Electromag const&   */ em,       //
            /*double mass         */ 1,
            /*Interpolator&*/ interpolator,
            /*GridLayout const&*/ layout, //
            no_op, no_op);
    }
}

#endif

#if PHARE_UPDATER_IMPL == 1

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    using namespace PHARE::core;

    auto static constexpr opts          = PHARE::SimOpts{dim, interp};
    bool static constexpr copy_particle = true;
    constexpr std::uint32_t cells       = 65;
    constexpr std::uint32_t n_parts     = 1e7;

    using PHARE_Types   = PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Electromag_t  = UsableElectromag<dim>;
    using Interpolator  = PHARE::core::Interpolator<dim, interp>;
    using ParticleArray = PHARE_Types::ParticleArray_t;
    using Particle_t    = ParticleArray::value_type;
    using Ions_t        = PHARE_Types::Ions_t;
    using Pusher        = IonUpdater1<Ions_t>::Pusher;
    using Boxing_t      = PHARE::core::UpdaterSelectionBoxing<GridLayout_t, ParticleArray>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    UsableIons<ParticleArray, interp> ions{layout};

    auto& domainParticles    = ions[0].domainParticles();
    domainParticles.vector() = std::vector<Particle_t>(n_parts, bench::particle<dim>());
    bench::disperse(domainParticles, 0, cells - 1, 13337);

    Boxing_t const boxing{layout, {grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts())}};
    Interpolator interpolator;
    Pusher pusher{.001, 1.0, layout};
    while (state.KeepRunningBatch(1))
        pusher.template move_interpolate_and_sort<true>(ions[0], em, boxing, interpolator);
}

#endif


BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
