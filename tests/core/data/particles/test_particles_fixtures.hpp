#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP


#include <core/data/particles/particle.hpp>
#include <core/data/ions/ion_population/particle_pack.hpp>

#include "test_particles.hpp"


#include <string>
#include <cassert>

namespace PHARE::core
{

template<typename ParticleArray_t>
struct UsableParticlesPopulation
{
    template<typename GridLayout_t>
    UsableParticlesPopulation(std::string const& _name, GridLayout_t const& layout)
        : name{_name}
        , domain_particles{make_particles<ParticleArray_t>(layout)}
    {
    }

    UsableParticlesPopulation(UsableParticlesPopulation const& that)
        : name{that.name}
        , domain_particles{that.domain_particles}
        , patch_ghost_particles{that.patch_ghost_particles}
        , level_ghost_particles{that.level_ghost_particles}
        , levelGhostParticlesOld{that.levelGhostParticlesOld}
        , levelGhostParticlesNew{that.levelGhostParticlesNew}
    {
    }

    auto& pack() { return particles_pack; }
    auto& pack() const { return particles_pack; }

    std::string name;
    ParticleArray_t domain_particles;
    ParticleArray_t patch_ghost_particles  = domain_particles;
    ParticleArray_t level_ghost_particles  = domain_particles;
    ParticleArray_t levelGhostParticlesOld = domain_particles;
    ParticleArray_t levelGhostParticlesNew = domain_particles;

    core::ParticlesPack<ParticleArray_t> particles_pack{name,
                                                        &domain_particles,
                                                        &patch_ghost_particles,
                                                        &level_ghost_particles,
                                                        &levelGhostParticlesOld,
                                                        &levelGhostParticlesNew};
};




} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP */
