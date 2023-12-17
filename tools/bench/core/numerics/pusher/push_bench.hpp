
#ifndef PHARE_BENCH_CORE_NUMERICS_PUSHER_PUSH_BENCH_HPP
#define PHARE_BENCH_CORE_NUMERICS_PUSHER_PUSH_BENCH_HPP


#include "tools/bench/core/bench.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"


namespace PHARE::core::bench
{



template<std::size_t dim>
void write_raw_unsorted_particles_to_file(std::size_t const n_parts = 1e7)
{
    constexpr std::uint32_t cells = 65;

    using PHARE_Types     = PHARE::core::PHARE_Types<dim, /*interp =*/1 /*not important here*/>;
    using GridLayout_t    = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using ParticleArray_t = typename PHARE_Types::ParticleArray_t;
    using Particle_t      = typename ParticleArray_t::value_type;

    GridLayout_t layout{cells};

    ParticleArray domainParticles{layout.AMRBox()};
    domainParticles.vector() = std::vector<Particle_t>(n_parts, core::bench::particle<dim>());
    core::bench::disperse(domainParticles, 0, cells - 1, 13337);

    std::stringstream ss;
    ss << "unsorted_particles_" << dim << ".raw";
    PHARE::core::bench::write_raw_to_file(domainParticles, ss.str());
}




} // namespace PHARE::core::bench

#endif /*PHARE_BENCH_CORE_NUMERICS_PUSHER_PUSH_BENCH_HPP*/
