#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP

#include <string>

namespace PHARE::core
{

template<typename ParticleArray_t>
struct UsableParticlesPopulation
{
    template<typename... Args>
    UsableParticlesPopulation(std::string const& _name, Args&&... args)
        : name{_name}
        , domain_particles{args...}
    {
    }

    auto& pack() { return particles_pack; }
    auto& pack() const { return particles_pack; }

    std::string name;
    ParticleArray_t domain_particles;
    ParticleArray_t patch_ghost_particles = domain_particles;
    ParticleArray_t level_ghost_particles = domain_particles;
    core::ParticlesPack<ParticleArray_t> particles_pack{name,
                                                        &domain_particles, //
                                                        &patch_ghost_particles,
                                                        &level_ghost_particles,
                                                        /*levelGhostParticlesOld=*/nullptr,
                                                        /*levelGhostParticlesNew=*/nullptr};
};


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP */
