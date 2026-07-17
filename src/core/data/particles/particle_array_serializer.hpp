#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERIALIZER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERIALIZER


#include "core/data/particles/serializing/particles_serializing.hpp"

namespace PHARE::core
{


template<typename Src>
void serialize_particles(std::string const& file_name, Src const& src)
{
    using Serializer = ParticlesSerializer<Src::layout_mode, Src::alloc_mode>;

    Serializer{}.operator()(file_name, src);
}

template<typename Dst, typename Src = Dst>
void deserialize_particles(std::string const& file_name, Dst& dst)
{
    using Deserializer = ParticlesDeserializer<Dst::layout_mode, Dst::alloc_mode>;

    Deserializer{}.template operator()<Dst, Src>(file_name, dst);
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SERIALIZER */
