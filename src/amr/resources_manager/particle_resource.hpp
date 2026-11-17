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
    template<typename ResourcesUser, typename GridLayout_t>
    struct UserParticleType
    {
        static constexpr auto dimension          = ResourcesUser::dimension;
        static constexpr std::size_t ghost_width = GridLayout_t::options.particle_ghost_width;

        using particle_array_type = ResourcesUser::particle_array_type;
        using variable_type       = ParticlesVariable<particle_array_type, ghost_width>;
        using patch_data_type     = ParticlesData<particle_array_type>;
    };

} // namespace amr
} // namespace PHARE


#endif // PARTICLE_RESOURCES_HPP
