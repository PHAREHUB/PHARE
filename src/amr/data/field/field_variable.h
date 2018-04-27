#ifndef PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_H
#define PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_H

#include <SAMRAI/hier/Variable.h>

#include <utility>

#include "data/grid/gridlayout.h"
#include "field_data_factory.h"

namespace PHARE
{
template<Layout layout, std::size_t dim, std::size_t interpOrder, typename FieldImpl,
         typename PhysicalQuantity = decltype(std::declval<FieldImpl>().physicalQuantity())>
class FieldVariable : public SAMRAI::hier::Variable
{
public:
    /** \brief Construct a new variable with an unique name, and a specific PhysicalQuantity
     *
     *  FieldVariable represent a data on a patch, it does not contain the data itself,
     *  after creation, one need to register it with a context : see registerVariableAndContext.
     */
    FieldVariable(std::string const& name, bool fineBoundaryRepresentsVariable,
                  PhysicalQuantity qty)
        : SAMRAI::hier::Variable(
              name,
              std::make_shared<FieldDataFactory<layout, dim, interpOrder, FieldImpl>>(
                  fineBoundaryRepresentsVariable, computeDataLivesOnPatchBorder_(qty), name, qty))
        , fineBoundaryRepresentsVariable_{fineBoundaryRepresentsVariable}
        , dataLivesOnPatchBorder_{computeDataLivesOnPatchBorder_(qty)}
    {
    }



    bool fineBoundaryRepresentsVariable() const final { return fineBoundaryRepresentsVariable_; }



    /** \brief Determines whether or not if data may lives on patch border
     *
     *  It will be true if in at least one direction, the data is primal
     */
    bool dataLivesOnPatchBorder() const final { return dataLivesOnPatchBorder_; }

private:
    bool const fineBoundaryRepresentsVariable_;
    bool const dataLivesOnPatchBorder_;



    bool computeDataLivesOnPatchBorder_(PhysicalQuantity qty)
    {
        auto const& centering = GridLayout<layout, dim>::centering(qty);


        for (auto const& qtyCentering : centering)
        {
            if (qtyCentering == QtyCentering::primal)
            {
                return true;
            }
        }
        return false;
    }
};



} // namespace PHARE
#endif
