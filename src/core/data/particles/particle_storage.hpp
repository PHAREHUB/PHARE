#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_STORAGE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_STORAGE_HPP


#include "core/utilities/multi_precision.hpp"


// Particle delta/velocity storage precision, see core/utilities/multi_precision.hpp

namespace PHARE::core
{
// Bytes == 8 is plain double
template<std::size_t Bytes>
struct particle_delta_storage
{
    using type = FixedPointUnit<Bytes>;
};
template<>
struct particle_delta_storage<8>
{
    using type = double;
};

template<std::size_t Bytes>
struct particle_velocity_storage
{
    using type = TruncatedDouble<Bytes>;
};
template<>
struct particle_velocity_storage<8>
{
    using type = double;
};


#ifndef PHARE_PARTICLE_DELTA_BYTES
#define PHARE_PARTICLE_DELTA_BYTES 6
#endif

#ifndef PHARE_PARTICLE_V_BYTES
#define PHARE_PARTICLE_V_BYTES 6
#endif

using ParticleDelta_t = particle_delta_storage<PHARE_PARTICLE_DELTA_BYTES>::type;
using ParticleV_t     = particle_velocity_storage<PHARE_PARTICLE_V_BYTES>::type;


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_STORAGE_HPP */
