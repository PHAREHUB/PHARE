#ifndef PHARE_PARTICLE_RESOURCES_H
#define PHARE_PARTICLE_RESOURCES_H



#include "core/data/particles/particle_array.h"
#include "amr/data/particles/particles_data.h"
#include "amr/data/particles/particles_variable.h"

namespace PHARE
{
namespace amr
{
    /** @brief tells SAMRAI which kind of variable, patchdata are used for a Field Resource
     * also says the type of the actual data buffer
     */
    template<class ResourcesUser>
    struct UserParticleType
    {
        static constexpr std::size_t dimension = ResourcesUser::dimension;

        using variable_type     = ParticlesVariable<dimension>;
        using patch_data_type   = ParticlesData<dimension>;
        using internal_type_ptr = typename ResourcesUser::particle_resource_type*;
    };

} // namespace amr
} // namespace PHARE


#endif // PARTICLE_RESOURCES_H
