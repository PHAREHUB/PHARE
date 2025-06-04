#ifndef PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP

#include <core/data/particles/particle.hpp>
#include <core/data/ions/ion_population/particle_pack.hpp>

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


template<std::size_t dim, typename Particle_t = Particle<dim>>
Particle_t particle(std::array<int, dim> const& icell)
{
    return {/*.weight = */ .001,
            /*.charge = */ 1,
            /*.iCell  = */ icell,
            /*.delta  = */ ConstArray<double, dim>(.51),
            /*.v      = */ {{.002002002002, .003003003003, .004004004004}}};
}

template<std::size_t dim>
Particle<dim> particle(int const icell = 15)
{
    return particle(ConstArray<int, dim>(icell));
}

template<typename Particles, typename Box>
void add_particles_in(Particles& particles, Box const& box, std::size_t const ppc)
{
    for (auto const& bix : box)
        for (std::size_t i = 0; i < ppc; ++i)
            particles.emplace_back(particle(*bix));
}


} // namespace PHARE::core


#endif /* PHARE_TEST_CORE_DATA_PARTICLES_TEST_PARTICLES_FIXTURES_HPP */
