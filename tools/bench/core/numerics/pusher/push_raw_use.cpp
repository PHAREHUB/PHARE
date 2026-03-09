#include "push_bench.hpp"
#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

// Does not include google benchmark as we want to see only PHARE operations/instructions/etc

template<std::size_t dim, std::size_t interp>
void push()
{
    auto static constexpr opts          = PHARE::SimOpts{dim, interp};
    bool static constexpr copy_particle = true;
    constexpr std::uint32_t cells       = 65;
    // constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types       = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t      = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Electromag_t      = PHARE::core::UsableElectromag<dim>;
    using Ions_t            = PHARE_Types::Ions_t;
    using ParticleArray     = Ions_t::particle_array_type;
    using IonUpdater        = typename PHARE::core::IonUpdater<Ions_t, Electromag_t, GridLayout_t>;
    using Boxing_t          = PHARE::core::UpdaterSelectionBoxing<IonUpdater, GridLayout_t>;

    using BorisPusher_t = PHARE::core::BorisPusher<dim, ParticleArray, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    Boxing_t boxing{layout, {layout.AMRBox()}};
    PHARE::core::UsableIons<ParticleArray, interp> ions{layout};

    std::stringstream ss;
    ss << "unsorted_particles_" << dim << ".raw";
    PHARE::core::bench::read_raw_from_file(ions[0].domainParticles(), ss.str());

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    Interpolator interpolator;
    pusher.template move<copy_particle>(ions[0], em, interpolator, boxing);
}

int main(int /*argc*/, char** /*argv*/)
{
    push<3, 3>();
}
