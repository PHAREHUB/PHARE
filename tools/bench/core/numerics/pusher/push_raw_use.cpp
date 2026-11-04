// PHARE_UPDATER_IMPL is set per target, see CMakeLists.txt
#include "push_bench.hpp"

#include "tests/core/data/ion_population/test_ion_population_fixtures.hpp"

// Does not include google benchmark as we want to see only PHARE operations/instructions/etc


namespace PHARE::core
{
template<std::size_t dim, std::size_t interp>
struct PushSetup
{
    auto static constexpr opts           = SimOpts{dim, interp};
    std::uint32_t static constexpr cells = 65;
    double static constexpr dt           = .001;
    double static constexpr charge       = 1;
    double static constexpr mass         = 1;

    using Types_t        = PHARE_Types<opts>;
    using GridLayout_t   = TestGridLayout<typename Types_t::Hybrid::GridLayout_t>;
    using Electromag_t   = UsableElectromag<dim>;
    using Interpolator_t = Interpolator<dim, interp>;
    using Ions_t         = Types_t::Hybrid::Ions_t;
    using ParticleArray  = Ions_t::particle_array_type;

    static auto raw_file() { return "unsorted_particles_" + std::to_string(dim) + ".raw"; }

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    Interpolator_t interpolator;
};
} // namespace PHARE::core



#if PHARE_UPDATER_IMPL == 0

template<std::size_t dim, std::size_t interp>
void push()
{
    using namespace PHARE::core;
    using Setup       = PushSetup<dim, interp>;
    using Particles_t = Setup::ParticleArray;
    using Pusher      = IonUpdater0<typename Setup::Ions_t>::Pusher;

    Setup s;
    Particles_t domainParticles{s.layout.AMRBox()};
    bench::read_raw_from_file(domainParticles, Setup::raw_file());

    Particles_t tmpDomain{s.layout.AMRBox()};
    tmpDomain.vector() = domainParticles.vector();

    auto rangeIn     = makeIndexRange(domainParticles);
    auto rangeOut    = makeIndexRange(tmpDomain);
    auto const no_op = [](auto& particleRange) { return particleRange; };

    Pusher const pusher{Setup::dt, Setup::charge, Setup::mass, s.layout};
    pusher.move(rangeIn, rangeOut, s.em, s.interpolator, no_op, no_op);
}

#endif


#if PHARE_UPDATER_IMPL == 1

template<std::size_t dim, std::size_t interp>
void push()
{
    using namespace PHARE::core;
    using Setup        = PushSetup<dim, interp>;
    using GridLayout_t = Setup::GridLayout_t;
    using Particles_t  = Setup::ParticleArray;
    using Pusher       = IonUpdater1<typename Setup::Ions_t>::Pusher;
    using Boxing_t     = UpdaterSelectionBoxing<GridLayout_t, Particles_t>;

    Setup s;
    UsableIons_t<Particles_t, interp> ions{s.layout};
    Boxing_t const boxing{s.layout,
                          {grow(s.layout.AMRBox(), GridLayout_t::options.particle_ghost_width)}};
    bench::read_raw_from_file(ions[0].domainParticles(), Setup::raw_file());

    Pusher const pusher{Setup::dt, Setup::charge, Setup::mass, s.layout};
    pusher.template move_interpolate_and_sort</*copy_particle=*/false>(ions[0], s.em, boxing,
                                                                       s.interpolator);
}

#endif



int main(int /*argc*/, char** /*argv*/)
{
    push<3, 3>();
}
