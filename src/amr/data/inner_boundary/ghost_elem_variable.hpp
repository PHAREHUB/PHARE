#ifndef PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_VARIABLE_HPP
#define PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_VARIABLE_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "ghost_elem_data_factory.hpp"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Variable.h>
#include <SAMRAI/tbox/Dimension.h>

#include <cstddef>
#include <string>

namespace PHARE::amr
{

template<std::size_t dim>
class GhostElemVariable : public SAMRAI::hier::Variable
{
public:
    explicit GhostElemVariable(std::string const& name)
        : SAMRAI::hier::Variable{name,
                                 std::make_shared<GhostElemDataFactory<dim>>(
                                     SAMRAI::hier::IntVector::getZero(
                                         SAMRAI::tbox::Dimension{dim}),
                                     name)}
    {
    }

    bool fineBoundaryRepresentsVariable() const final { return true; }
    bool dataLivesOnPatchBorder() const final { return false; }
};

} // namespace PHARE::amr

#endif // PHARE_AMR_DATA_INNER_BOUNDARY_GHOST_ELEM_VARIABLE_HPP
