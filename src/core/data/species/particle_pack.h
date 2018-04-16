#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_PACK_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_PACK_H

#include "data/particles/particle_array.h"

namespace PHARE
{
//! ParticlePack is used to pack three ParticleArray together
/**
 * For conveniency, particles are often stored in arrays specific to their
 * location. Ghost particles are particles coming from neighbor patches of
 * the same level, incoming particles are particles coming from coarse-to-fine
 * boundaries, and interior particles are those within the patch physical domain.
 * ParticlesPack conveniently store the three arrays together.
 */
template<std::size_t dim>
struct ParticlesPack
{
    ParticleArray<dim> interior;
    ParticleArray<dim> ghost;
    ParticleArray<dim> incoming;
};

} // namespace PHARE

#endif // PARTICLE_PACK_H
