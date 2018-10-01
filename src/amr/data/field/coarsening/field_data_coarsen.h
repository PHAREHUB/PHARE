#ifndef PHARE_FIELD_DATA_COARSEN_H
#define PHARE_FIELD_DATA_COARSEN_H

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/IntVector.h>

#include "data/field/field_data.h"
#include "data/field/field_geometry.h"
#include "field_coarsen.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"

namespace PHARE
{
//
//
template<typename GridLayoutT, typename FieldT,
         typename PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity())>
class FieldDataCoarsen : public SAMRAI::hier::CoarsenOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    static constexpr std::size_t maxRafinement{10};
    FieldDataCoarsen()
        : SAMRAI::hier::CoarsenOperator("FieldDataCoarsenOperator")
    {
    }

    FieldDataCoarsen(FieldDataCoarsen const&) = delete;
    FieldDataCoarsen(FieldDataCoarsen&&)      = delete;
    FieldDataCoarsen& operator=(FieldDataCoarsen const&) = delete;
    FieldDataCoarsen&& operator=(FieldDataCoarsen&&) = delete;


    virtual ~FieldDataCoarsen() = default;




    /** @brief return the priority of the operator
     *  this return 0, meaning that this operator
     * have the most priority
     */
    int getOperatorPriority() const override { return 0; }




    /** @brief Return the stencil width associated with the coarsening operator.
     *
     *  The SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghostCell data surrounding the interior to satisfy the stencil
     * width requirements for each coarsening operator.
     *
     * In our case, we allow a RF up to 10, so having 5 ghost width is sufficient
     *
     */
    SAMRAI::hier::IntVector getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector{dim, maxRafinement / 2};
    }




    /** @brief given a coarseBox, coarse data from the fine patch on the intersection of this box
     * and the box of the destination (the box of the coarse patch).
     *
     * This method will extract fieldData from the two patches, and then
     * get the Field and GridLayout encapsulated into the fieldData.
     * With the help of FieldGeometry, transform the coarseBox to the correct index.
     * After that we can now create FieldCoarsen with the indexAndWeight implementation selected.
     * Finnaly loop over the indexes in the box, and apply the coarsening defined in FieldCoarsen
     * operator
     *
     */
    void coarsen(SAMRAI::hier::Patch& destinationPatch, SAMRAI::hier::Patch const& sourcePatch,
                 int const destinationComponent, int const sourceComponent,
                 SAMRAI::hier::Box const& coarseBox,
                 SAMRAI::hier::IntVector const& ratio) const override
    {
        auto destinationFieldData = std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldT>>(
            destinationPatch.getPatchData(destinationComponent));
        auto const sourceFieldData
            = std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldT> const>(
                sourcePatch.getPatchData(sourceComponent));

        if (!destinationFieldData || !sourceFieldData)
        {
            throw std::runtime_error("cannot to FieldData");
        }

        // We get layout from the fieldData
        auto const& destinationLayout = destinationFieldData->gridLayout;
        auto const& sourceLayout      = sourceFieldData->gridLayout;

        // We get field from the fieldData
        auto& destinationField  = destinationFieldData->field;
        auto const& sourceField = sourceFieldData->field;

        // we assume that quantity are the same
        // note that an assertion will be raised
        // in coarseIt operator
        auto const& qty = destinationField.physicalQuantity();


        bool const withGhost{true};

        // We get different boxes : destination , source, restrictBoxes
        // and transform them in the correct indexing.
        auto destinationBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            destinationPatch.getBox(), qty, destinationLayout, withGhost);

        auto sourceBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            sourcePatch.getBox(), qty, sourceLayout, withGhost);

        auto coarseLayout = FieldGeometry<GridLayoutT, PhysicalQuantity>::layoutFromBox(
            coarseBox, destinationLayout);

        auto coarseFieldBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            coarseBox, qty, coarseLayout, !withGhost);

        // finnaly we compute the intersection
        auto intersectionBox = destinationBox * coarseFieldBox;




        // We can now create the coarsening operator
        CoarsenField<dimension> coarsenIt{destinationLayout.centering(qty), sourceBox,
                                          destinationBox, ratio};

        // now we can loop over the intersection box

        Point<int, dimension> startIndex;
        Point<int, dimension> endIndex;

        startIndex[dirX] = intersectionBox.lower(dirX);
        endIndex[dirX]   = intersectionBox.upper(dirX);

        if constexpr (dimension > 1)
        {
            startIndex[dirY] = intersectionBox.lower[dirY];
            endIndex[dirY]   = intersectionBox.upper[dirY];
        }
        if constexpr (dimension > 2)
        {
            startIndex[dirZ] = intersectionBox.lower[dirZ];
            endIndex[dirZ]   = intersectionBox.upper[dirZ];
        }

        if constexpr (dimension == 1)
        {
            for (int ix = startIndex[dirX]; ix <= endIndex[dirX]; ++ix)
            {
                coarsenIt(sourceField, destinationField, {{ix}});
            }
        }
        else if constexpr (dimension == 2)
        {
            for (int ix = startIndex[dirX]; ix <= endIndex[dirX]; ++ix)
            {
                for (int iy = startIndex[dirY]; iy <= endIndex[dirY]; ++iy)
                {
                    coarsenIt(sourceField, destinationField, {{ix, iy}});
                }
            }
        }
        else if constexpr (dimension == 3)
        {
            for (int ix = startIndex[dirX]; ix <= endIndex[dirX]; ++ix)
            {
                for (int iy = startIndex[dirY]; iy <= endIndex[dirY]; ++iy)
                {
                    for (int iz = startIndex[dirZ]; iz <= endIndex[dirZ]; ++iz)

                    {
                        coarsenIt(sourceField, destinationField, {{ix, iy, iz}});
                    }
                }
            }
        }
    }
};
} // namespace PHARE


#endif
