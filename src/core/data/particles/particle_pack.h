#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_PACK_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_PACK_H

#include "particle_array.h"

//! ParticlePack is used to pack three ParticleArray together
/**
 *
 */
class ParticlesPack
{
public:
    ParticleArray interior;
    ParticleArray ghost;
    ParticleArray incoming;
};


#endif // PARTICLE_PACK_H
