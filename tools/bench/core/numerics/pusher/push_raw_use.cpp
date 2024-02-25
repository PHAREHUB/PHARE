#include "push_bench.hpp"

// Does not include google benchmark as we want to see only PHARE operations/instructions/etc

template<std::size_t dim, std::size_t interp>
void push()
{
    constexpr std::uint32_t cells = 65;
    // constexpr std::uint32_t n_parts = 1e7;

    using PHARE_Types       = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t      = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Electromag_t      = PHARE::core::bench::Electromag<GridLayout_t>;
    using Ions_t            = typename PHARE_Types::Ions_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using ParticleRange     = PHARE::core::IndexRange<ParticleArray>;
    using BorisPusher_t = PHARE::core::BorisPusher<dim, ParticleRange, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};

    ParticleArray domainParticles{layout.AMRBox()};

    std::stringstream ss;
    ss << "unsorted_particles_" << dim << ".raw";
    PHARE::core::bench::read_raw_from_file(domainParticles, ss.str());

    // std::sort(domainParticles);

    ParticleArray tmpDomain{layout.AMRBox()};
    tmpDomain.vector() = domainParticles.vector();

    auto rangeIn  = PHARE::core::makeIndexRange(domainParticles);
    auto rangeOut = PHARE::core::makeIndexRange(tmpDomain);

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    Interpolator interpolator;
    auto const no_op = [](auto& particleRange) { return particleRange; };

    pusher.move(
        /*ParticleRange const&*/ rangeIn,  //
        /*ParticleRange&      */ rangeOut, //
        /*Electromag const&   */ em,       //
        /*double mass         */ 1,
        /*Interpolator&*/ interpolator,
        /*GridLayout const&*/ layout, //
        no_op, no_op);
}

int main(int /*argc*/, char** /*argv*/)
{
    push<3, 3>();
}
