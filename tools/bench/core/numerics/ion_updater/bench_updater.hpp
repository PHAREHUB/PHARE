
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
        return IonUpdaterProxy<IonUpdater0<Ions>>();
    if constexpr (impl == 1)
        return IonUpdaterProxy<IonUpdater1<Ions>>();
    if constexpr (impl == 2)
        return IonUpdaterProxy<IonUpdater2<Ions>>();
}

template<std::size_t dim_, std::size_t interp_, std::size_t impl_>
struct Params
{
    static constexpr auto dim    = dim_;
    static constexpr auto interp = interp_;
    static constexpr auto impl   = impl_;
};

template<typename Params>
struct UpdaterBencher : public benchmark::Fixture
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
    using Ions          = UsableIons<ParticleArray, interp>;
    // using IonUpdater    = IonUpdaterProxy<IonUpdater2<Ions>>;
    using Boxing_t = UpdaterSelectionBoxing<GridLayout_t, ParticleArray>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    Ions ions{layout};
    Boxing_t const boxing{layout, {grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts())}};
    ParticleArray particles_copy;

    void SetUp(::benchmark::State& state)
    {
        auto& patch_particles = ions.populations[0].particles;
        patch_particles.domain_particles.vector()
            = std::vector<Particle_t>(n_parts, bench::particle<dim>());
        bench::disperse(patch_particles.domain_particles, 0, cells - 1);
        std::sort(patch_particles.domain_particles);
        particles_copy = patch_particles.domain_particles; // tmp storage between update modes
    }

    void TearDown(::benchmark::State& state) {}

    void operator()(benchmark::State& state)
    {
        initializer::PHAREDict dict;
        dict["pusher"]["name"] = std::string{"modified_boris"};
        auto ionUpdater_       = get_updater_for<impl>(ions);

        auto& patch_particles     = ions.populations[0].particles;
        double const current_time = 1.0;
        double const new_time     = 1.005;
        auto const dt             = new_time - current_time;
        while (state.KeepRunningBatch(1)) // while (state.KeepRunning())
        {
            ionUpdater_.updatePopulations(UpdaterMode::domain_only, ions, em, boxing, dt);
            ionUpdater_.updateIons(ions);

            patch_particles.domain_particles = particles_copy;
            auto& pack                       = std::get<4>(
                ions.getRunTimeResourcesViewList()[0].getCompileTimeResourcesViewList());
            pack.setBuffer(&patch_particles.pack());


            ionUpdater_.updatePopulations(UpdaterMode::all, ions, em, boxing, dt);
            ionUpdater_.updateIons(ions);
        }
    }
};

} // namespace PHARE::core
