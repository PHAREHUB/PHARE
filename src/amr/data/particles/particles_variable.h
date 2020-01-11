#ifndef PHARE_PARTICLES_VARIABLE_H
#define PHARE_PARTICLES_VARIABLE_H

#include "particles_data_factory.h"
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/tbox/Dimension.h>

namespace PHARE
{
namespace amr
{
    /**
     * @brief The ParticlesVariable class
     */
    template<std::size_t dim, std::size_t interp>
    class ParticlesVariable : public SAMRAI::hier::Variable
    {
    public:
        ParticlesVariable(std::string const& name, bool fineBoundaryRepresentsVariable = false,
                          SAMRAI::hier::IntVector ghost
                          = SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dim},
                                                    ghostWidthForParticles<interp>()})
            : SAMRAI::hier::Variable{name, std::make_shared<ParticlesDataFactory<dim>>(
                                               ghost, fineBoundaryRepresentsVariable)}
            , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        {
        }




        virtual bool fineBoundaryRepresentsVariable() const final
        {
            return fineBoundaryRepresentsVariable_;
        }



        virtual bool dataLivesOnPatchBorder() const final { return true; }

    private:
        bool fineBoundaryRepresentsVariable_;
    };

} // namespace amr

} // namespace PHARE
#endif
