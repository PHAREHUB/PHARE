
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "phare_core.hpp"
#include "phare_simulator_options.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

#include "tools/bench/core/bench.hpp"

#include "benchmark/benchmark.h"

namespace PHARE::core
{

template<std::uint16_t impl, typename Ions>
auto get_updater_for(Ions const&)
{
    if constexpr (impl == 0)
        return IonUpdater0<Ions>();
    if constexpr (impl == 1)
        return IonUpdater1<Ions>();
}

template<std::size_t dim_, std::size_t interp_, std::size_t impl_>
struct Params
{
    static constexpr auto dim    = dim_;
    static constexpr auto interp = interp_;
    static constexpr auto impl   = impl_;
};


template<typename Params>
struct BenchSetup
{
    static constexpr auto dim              = Params::dim;
    static constexpr auto interp           = Params::interp;
    static constexpr auto impl             = Params::impl;
    static constexpr std::uint32_t cells   = 30;
    static constexpr std::uint32_t n_parts = 1e7;
    static constexpr auto opts             = PHARE::SimOpts{dim, interp};

    using PHARE_Types   = core::PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Electromag_t  = UsableElectromag<dim>;
    using ParticleArray = PHARE_Types::ParticleArray_t;
    using Particle_t    = ParticleArray::value_type;
    using Ions          = UsableIons_t<ParticleArray, interp>;
    using Boxing_t      = UpdaterSelectionBoxing<GridLayout_t, ParticleArray>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    Ions ions{layout};
    Boxing_t const boxing{layout, {grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts())}};
    ParticleArray particles_copy;

    void SetUp(::benchmark::State&)
    {
        auto& patch_particles = ions.populations[0].particles;
        patch_particles.domain_particles.vector()
            = std::vector<Particle_t>(n_parts, bench::particle<dim>());
        bench::disperse(patch_particles.domain_particles, 0, cells - 1, 133337);
        std::sort(patch_particles.domain_particles);
        patch_particles.domain_particles.map_particles();
        particles_copy = patch_particles.domain_particles;
    }
};


// Benchmark UpdaterMode::copy only.
// Domain particles are not mutated, so no restoration between iterations.
template<typename Params>
struct UpdaterBencherCopy : public benchmark::Fixture, public BenchSetup<Params>
{
    using Setup = BenchSetup<Params>;

    void SetUp(::benchmark::State& state) override { Setup::SetUp(state); }
    void TearDown(::benchmark::State&) override {}

    void operator()(benchmark::State& state)
    {
        auto ionUpdater_ = get_updater_for<Setup::impl>(Setup::ions);
        auto const dt    = 0.005;
        while (state.KeepRunningBatch(1))
        {
            ionUpdater_.updatePopulations(UpdaterMode::copy, Setup::ions, Setup::em, Setup::boxing,
                                          dt);
            Setup::ions.update();
        }
    }
};


// Benchmark UpdaterMode::ref only.
// Domain particles move each iteration, so restore before each run.
template<typename Params>
struct UpdaterBencherRef : public benchmark::Fixture, public BenchSetup<Params>
{
    using Setup = BenchSetup<Params>;

    void SetUp(::benchmark::State& state) override { Setup::SetUp(state); }
    void TearDown(::benchmark::State&) override {}

    void operator()(benchmark::State& state)
    {
        auto ionUpdater_ = get_updater_for<Setup::impl>(Setup::ions);
        auto const dt    = 0.005;
        while (state.KeepRunningBatch(1))
        {
            state.PauseTiming();
            auto& patch_particles            = Setup::ions.populations[0].particles;
            patch_particles.domain_particles = Setup::particles_copy;
            auto& pack                       = std::get<4>(
                Setup::ions.getRunTimeResourcesViewList()[0].getCompileTimeResourcesViewList());
            pack.setBuffer(&patch_particles.pack());
            state.ResumeTiming();

            ionUpdater_.updatePopulations(UpdaterMode::ref, Setup::ions, Setup::em, Setup::boxing,
                                          dt);
            Setup::ions.update();
        }
    }
};

} // namespace PHARE::core
