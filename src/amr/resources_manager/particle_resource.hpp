#ifndef PHARE_PARTICLE_RESOURCES_HPP
#define PHARE_PARTICLE_RESOURCES_HPP


#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_variable.hpp"

namespace PHARE
{
namespace amr
{
    /** @brief tells SAMRAI which kind of variable, patchdata are used for a Particle resource
     */
    template<typename ParticleArray, std::size_t interp>
    struct ParticleViewInfo
    {
        static constexpr auto dimension    = ParticleArray::dimension;
        static constexpr auto interp_order = interp;

        // TODORM this class could be templated by the ParticleArray type and rather
        // hard code the link to the ParticlePack here
        using particle_array_type = ParticleArray;
        using variable_type       = ParticlesVariable<particle_array_type, interp_order>;
        using patch_data_type     = ParticlesData<particle_array_type>;
        using view_type           = patch_data_type::view_type;
    };

} // namespace amr
} // namespace PHARE


#endif // PARTICLE_RESOURCES_HPP
