#ifndef PHARE_PARTICLES_VARIABLE_H
#define PHARE_PARTICLES_VARIABLE_H

#include <SAMRAI/hier/Variable.h>

#include "particles_data_factory.h"

namespace PHARE
{
template<std::size_t dim>
class ParticlesVariable : public SAMRAI::hier::Variable
{
public:
    ParticlesVariable(std::string const& name, bool fineBoundaryRepresentsVariable,
                      SAMRAI::hier::IntVector ghost)
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



} // namespace PHARE
#endif
