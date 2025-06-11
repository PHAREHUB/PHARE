#ifndef PHARE_PARTICLES_VARIABLE_HPP
#define PHARE_PARTICLES_VARIABLE_HPP

#include "core/def/phare_mpi.hpp"        // IWYU pragma: keep
#include "core/data/grid/gridlayout.hpp" // particle ghost width

#include "particles_data_factory.hpp"

#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/tbox/Dimension.h>

namespace PHARE
{
namespace amr
{
    /**
     * @brief The ParticlesVariable class
     */
    template<typename ParticleArray, std::size_t interp>
    class ParticlesVariable : public SAMRAI::hier::Variable
    {
        static constexpr auto dim = ParticleArray::dimension;

    public:
        ParticlesVariable(std::string const& name, bool fineBoundaryRepresentsVariable = false,
                          SAMRAI::hier::IntVector ghost
                          = SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dim},
                                                    core::ghostWidthForParticles<interp>()})
            : SAMRAI::hier::Variable{name, std::make_shared<ParticlesDataFactory<ParticleArray>>(
                                               ghost, fineBoundaryRepresentsVariable, name)}
            , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        {
        }




        bool fineBoundaryRepresentsVariable() const final
        {
            return fineBoundaryRepresentsVariable_;
        }



        bool dataLivesOnPatchBorder() const final { return true; }

    private:
        bool fineBoundaryRepresentsVariable_;
    };

} // namespace amr

} // namespace PHARE
#endif
