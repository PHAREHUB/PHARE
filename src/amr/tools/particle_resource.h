#ifndef PHARE_PARTICLE_RESOURCES_H
#define PHARE_PARTICLE_RESOURCES_H



#include "data/particles/particle_array.h"
#include "data/particles/particles_data.h"
#include "data/particles/particles_variable.h"

namespace PHARE
{
/** @brief tells SAMRAI which kind of variable, patchdata are used for a Field Resource
 * also says the type of the actual data buffer
 */
template<typename ResourcesUser>
struct UserParticleType
{
    using patch_data_type   = ParticlesData<ResourcesUser::dimension>;
    using variable_type     = ParticlesVariable<ResourcesUser::dimension>;
    using internal_type_ptr = typename ResourcesUser::particle_array_type*;
};


} // namespace PHARE


#endif // PARTICLE_RESOURCES_H
