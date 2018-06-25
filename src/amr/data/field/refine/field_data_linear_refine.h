#ifndef PHARE_FIELD_DATA_LINEAR_REFINE_H
#define PHARE_FIELD_DATA_LINEAR_REFINE_H

#include "data/field/field_data.h"
#include "data/field/field_geometry.h"
#include "field_linear_refine.h"

#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/tbox/Dimension.h>

#include <string>


namespace PHARE
{
template<typename GridLayoutImpl, typename FieldT,
         typename PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity())>
class FieldDataLinearRefine : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutImpl::dimension;

    FieldDataLinearRefine()
        : SAMRAI::hier::RefineOperator{"FieldDataLinearRefineOperator"}
    {
    }

    virtual ~FieldDataLinearRefine() = default;

    /** This implementation have the top priority for refine operation
     *
     */
    int getOperatorPriority() const override { return 0; }

    /**
     * @brief This operator need to have at least 1 ghost cell to work properly
     *
     */
    SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getOne(dim);
    }

    /**
     * @brief Given a  set of box on a  fine patch, compute the interpolation from
     * a coarser patch that is underneath the fine box.
     * Since we get our boxes from a FieldOverlap, we know that there are in correct
     * Field Indexes
     *
     */
    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationComponent, int const sourceComponent,
                SAMRAI::hier::BoxOverlap const& destinationOverlap,
                SAMRAI::hier::IntVector const& ratio) const override
    {
        // cast the BoxOverlap to FieldOverlap ones;
        auto const& destinationFieldOverlap
            = dynamic_cast<FieldOverlap<dimension> const&>(destinationOverlap);

        auto const& destinationBoxes = destinationFieldOverlap.getDestinationBoxContainer();

        auto destinationFieldData = std::dynamic_pointer_cast<FieldData<GridLayoutImpl, FieldT>>(
            destination.getPatchData(destinationComponent));

        auto const sourceFieldData = std::dynamic_pointer_cast<FieldData<GridLayoutImpl, FieldT>>(
            source.getPatchData(sourceComponent));

        if (!destinationFieldData || !sourceFieldData)
        {
            throw std::runtime_error("Cannot cast to FieldData");
        }

        // We get layout from the fieldData
        auto const& destinationLayout = destinationFieldData->gridLayout;
        auto const& sourceLayout      = sourceFieldData->gridLayout;

        // We get field from fieldData
        auto& destinationField  = destinationFieldData->field;
        auto const& sourceField = sourceFieldData->field;

        // We assume that quantity are the same
        // note that an assertion will be raised
        // in refineIt operator
        auto const& qty = destinationField.physicalQuantity();

        bool const withGhost{true};

        auto destinationFieldBox = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
            destination.getBox(), qty, destinationLayout, withGhost);

        auto sourceFieldBox = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
            source.getBox(), qty, sourceLayout, withGhost);

        FieldLinearRefine<dimension> refineIt{destinationLayout.centering(qty), destinationFieldBox,
                                              sourceFieldBox, ratio};

        for (auto const& box : destinationBoxes)
        {
            // we compute the intersection with the destination,
            // and then we apply the refine operation on each fine
            // index.
            auto intersectionBox = destinationFieldBox * box;

            if constexpr (dimension == 1)
            {
                int iStartX = intersectionBox.lower(dirX);
                int iEndX   = intersectionBox.upper(dirX);

                for (int ix = iStartX; ix <= iEndX; ++ix)
                {
                    refineIt(sourceField, destinationField, {{ix}});
                }
            }
            else if constexpr (dimension == 2)
            {
                int iStartX = intersectionBox.lower(dirX);
                int iStartY = intersectionBox.lower(dirY);

                int iEndX = intersectionBox.upper(dirX);
                int iEndY = intersectionBox.upper(dirY);

                for (int ix = iStartX; ix <= iEndX; ++ix)
                {
                    for (int iy = iStartY; iy <= iEndY; ++iy)
                    {
                        refineIt(sourceField, destinationField, {{ix, iy}});
                    }
                }
            }
            else if constexpr (dimension == 3)
            {
                int iStartX = intersectionBox.lower(dirX);
                int iStartY = intersectionBox.lower(dirY);
                int iStartZ = intersectionBox.lower(dirZ);

                int iEndX = intersectionBox.upper(dirX);
                int iEndY = intersectionBox.upper(dirY);
                int iEndZ = intersectionBox.upper(dirZ);

                for (int ix = iStartX; ix <= iEndX; ++ix)
                {
                    for (int iy = iStartY; iy <= iEndY; ++iy)
                    {
                        for (int iz = iStartZ; iz <= iEndZ; ++iz)
                        {
                            refineIt(sourceField, destinationField, {{ix, iy, iz}});
                        }
                    }
                }
            }
        }
    }

private:
};

} // namespace PHARE



#endif
