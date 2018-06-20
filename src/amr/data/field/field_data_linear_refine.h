#ifndef PHARE_FIELD_DATA_LINEAR_REFINE_H
#define PHARE_FIELD_DATA_LINEAR_REFINE_H

#include "data/refine/field_linear_refine.h"
#include "field_data.h"
#include "field_geometry.h"

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
                int xLower = intersectionBox.lower(dirX);
                int xUpper = intersectionBox.upper(dirX);

                for (int iStartX = xLower; iStartX <= xUpper; ++iStartX)
                {
                    refineIt(sourceField, destinationField, {{iStartX}});
                }
            }
            else if constexpr (dimension == 2)
            {
                int xLower = intersectionBox.lower(dirX);
                int yLower = intersectionBox.lower(dirY);

                int xUpper = intersectionBox.upper(dirX);
                int yUpper = intersectionBox.upper(dirY);

                for (int iStartX = xLower; iStartX <= xUpper; ++iStartX)
                {
                    for (int iStartY = yLower; iStartY <= yUpper; ++iStartY)
                    {
                        refineIt(sourceField, destinationField, {{iStartX, iStartY}});
                    }
                }
            }
            else if constexpr (dimension == 3)
            {
                int xLower = intersectionBox.lower(dirX);
                int yLower = intersectionBox.lower(dirY);
                int zLower = intersectionBox.lower(dirZ);

                int xUpper = intersectionBox.upper(dirX);
                int yUpper = intersectionBox.upper(dirY);
                int zUpper = intersectionBox.upper(dirZ);

                for (int iStartX = xLower; iStartX <= xUpper; ++iStartX)
                {
                    for (int iStartY = yLower; iStartY <= yUpper; ++iStartY)
                    {
                        for (int iStartZ = zLower; iStartZ <= zUpper; ++iStartZ)
                        {
                            refineIt(sourceField, destinationField, {{iStartX, iStartY, iStartZ}});
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
