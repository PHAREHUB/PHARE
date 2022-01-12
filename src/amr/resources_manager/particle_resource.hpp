#ifndef PHARE_PARTICLE_RESOURCES_HPP
#define PHARE_PARTICLE_RESOURCES_HPP


#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_variable.hpp"

namespace PHARE
{
namespace amr
{
    /** @brief tells SAMRAI which kind of variable, patchdata are used for a Field Resource
     * also says the type of the actual data buffer
     */
    template<typename ResourcesUser, std::size_t interp>
    struct UserParticleType
    {
        static constexpr auto dimension    = ResourcesUser::dimension;
        static constexpr auto interp_order = interp;

        using particle_array_type = typename ResourcesUser::particle_array_type;
        using variable_type       = ParticlesVariable<particle_array_type, interp_order>;
        using patch_data_type     = ParticlesData<particle_array_type>;
        using internal_type_ptr   = typename ResourcesUser::particle_resource_type*;
    };

} // namespace amr
} // namespace PHARE


#endif // PARTICLE_RESOURCES_HPP
