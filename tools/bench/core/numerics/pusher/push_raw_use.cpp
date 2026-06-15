#include "core/numerics/ion_updater/ion_updater/ion_updater_impl0.hpp"
#define PHARE_UPDATER_IMPL 0

#include "push_bench.hpp"

#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

// Does not include google benchmark as we want to see only PHARE operations/instructions/etc



#if PHARE_UPDATER_IMPL == 0

template<std::size_t dim, std::size_t interp>
void push()
{
    using namespace PHARE::core;

    auto static constexpr opts    = PHARE::SimOpts{dim, interp};
    constexpr std::uint32_t cells = 65;
    // constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types   = PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator  = Interpolator<dim, interp>;
    using Electromag_t  = UsableElectromag<dim>;
    using Ions_t        = PHARE_Types::Ions_t;
    using ParticleArray = Ions_t::particle_array_type;
    using Pusher        = IonUpdater0<Ions_t>::Pusher;
    // using ParticleRange     = IndexRange<ParticleArray>;
    // using BoundaryCondition = BoundaryCondition<dim, interp>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};

    ParticleArray domainParticles{layout.AMRBox()};

    std::stringstream ss;
    ss << "unsorted_particles_" << dim << ".raw";
    bench::read_raw_from_file(domainParticles, ss.str());

    // std::sort(domainParticles);

    ParticleArray tmpDomain{layout.AMRBox()};
    tmpDomain.vector() = domainParticles.vector();

    auto rangeIn  = makeIndexRange(domainParticles);
    auto rangeOut = makeIndexRange(tmpDomain);

    Interpolator interpolator;
    auto const no_op = [](auto& particleRange) { return particleRange; };

    Pusher const pusher{/*dt=*/.001, /*mass=*/1, layout};
    pusher.move(
        /*ParticleRange const&*/ rangeIn,  //
        /*ParticleRange&      */ rangeOut, //
        /*Electromag const&   */ em,       //
        /*Interpolator&*/ interpolator, no_op, no_op);
}

#endif


#if PHARE_UPDATER_IMPL == 1

template<std::size_t dim, std::size_t interp>
void push()
{
    using namespace PHARE::core;

    auto static constexpr opts          = PHARE::SimOpts{dim, interp};
    bool static constexpr copy_particle = true;
    constexpr std::uint32_t cells       = 65;
    // constexpr std::uint32_t n_parts = 1e7;
    using PHARE_Types   = PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Electromag_t  = UsableElectromag<dim>;
    using Interpolator  = Interpolator<dim, interp>;
    using Ions_t        = PHARE_Types::Ions_t;
    using ParticleArray = Ions_t::particle_array_type;
    using Pusher        = IonUpdater1<Ions_t>::Pusher;
    using UsableIons_t  = UsableIons<ParticleArray, interp>;
    using Boxing_t      = PHARE::core::UpdaterSelectionBoxing<GridLayout_t, ParticleArray>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    UsableIons_t ions{layout};
    Boxing_t const boxing{layout, {grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts())}};

    std::stringstream ss;
    ss << "unsorted_particles_" << dim << ".raw";
    bench::read_raw_from_file(ions[0].domainParticles(), ss.str());

    Interpolator interpolator;
    Pusher const pusher{/*dt=*/.001, /*mass=*/1, layout};
    pusher.move(ions[0].domainParticles(), em, boxing, interpolator);
}

#endif



int main(int /*argc*/, char** /*argv*/)
{
    push<3, 3>();
}
