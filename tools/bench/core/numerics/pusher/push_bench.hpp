
#ifndef PHARE_BENCH_CORE_NUMERICS_PUSHER_PUSH_BENCH_HPP
#define PHARE_BENCH_CORE_NUMERICS_PUSHER_PUSH_BENCH_HPP


#include "tools/bench/core/bench.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"


namespace PHARE::core::bench
{



template<std::size_t dim>
void write_raw_unsorted_particles_to_file(std::size_t const n_parts = 1e6)
{
    constexpr std::uint32_t cells = 15;

    using PHARE_Types     = core::PHARE_Types<SimOpts{dim, /*interp =*/1 /*not important here*/}>;
    using GridLayout_t    = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using ParticleArray_t = PHARE_Types::ParticleArray_t;

    GridLayout_t layout{cells};

    ParticleArray_t domainParticles{layout.AMRBox()};
    add_particles_in(domainParticles, layout.AMRBox(), n_parts / layout.AMRBox().size());
    shuffle(domainParticles, 13337);

    std::stringstream ss;
    ss << "unsorted_particles_" << dim << ".raw";
    PHARE::core::write_raw_to_file(domainParticles, ss.str());
}




} // namespace PHARE::core::bench

#endif /*PHARE_BENCH_CORE_NUMERICS_PUSHER_PUSH_BENCH_HPP*/
