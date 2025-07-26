#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP


#include "core/def/phare_mpi.hpp"

#include "core/def.hpp"
#include "amr/data/field/field_data.hpp"

#include "field_linear_refine.hpp"

#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/RefineOperator.h>


#include <cstddef>


namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;

template<typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
class FieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = typename GridLayoutT::implT;
    using PhysicalQuantity                 = typename FieldT::physical_quantity_type;
    using FieldDataT                       = FieldData<GridLayoutT, FieldT>;

    FieldRefineOperator(bool node_only = false)
        : SAMRAI::hier::RefineOperator{"FieldRefineOperator"}

    {
    }

    virtual ~FieldRefineOperator() = default;

    /** This implementation have the top priority for refine operation
     *
     */
    NO_DISCARD int getOperatorPriority() const override { return 0; }

    /**
     * @brief This operator needs to have at least 1 ghost cell to work properly
     *
     */
    NO_DISCARD SAMRAI::hier::IntVector
    getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getOne(dim);
    }




    /**
     * @brief Given a set of box on a fine patch, compute the interpolation from
     * a coarser patch that is underneath the fine box.
     * Since we get our boxes from a FieldOverlap, we know that they are in correct
     * Field Indexes
     *
     */
    void refine(SAMRAI::hier::Patch& destination, SAMRAI::hier::Patch const& source,
                int const destinationId, int const sourceId,
                SAMRAI::hier::BoxOverlap const& destinationOverlap,
                SAMRAI::hier::IntVector const& ratio) const override
    {
        using FieldGeometry = typename FieldDataT::Geometry;

        auto const& destinationFieldOverlap = dynamic_cast<FieldOverlap const&>(destinationOverlap);

        auto const& overlapBoxes = destinationFieldOverlap.getDestinationBoxContainer();

        auto& destinationField  = FieldDataT::getField(destination, destinationId);
        auto const& destLayout  = FieldDataT::getLayout(destination, destinationId);
        auto const& sourceField = FieldDataT::getField(source, sourceId);
        auto const& srcLayout   = FieldDataT::getLayout(source, sourceId);


        // We assume that quantity are all the same.
        // Note that an assertion will be raised
        // in refineIt operator
        auto const& qty = destinationField.physicalQuantity();


        auto const destData = destination.getPatchData(destinationId);
        auto const srcData  = source.getPatchData(sourceId);

        auto const destFieldBox
            = FieldGeometry::toFieldBox(destData->getGhostBox(), qty, destLayout);
        auto const sourceFieldBox
            = FieldGeometry::toFieldBox(srcData->getGhostBox(), qty, srcLayout);


        FieldRefinerPolicy refiner{destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio};


        for (auto const& box : overlapBoxes)
        {
            // we compute the intersection with the destination,
            // and then we apply the refine operation on each fine
            // index.
            auto intersectionBox = destFieldBox * box;


            std::cout << "debug refine operator: "
                      << "destinationFieldBox: " << destFieldBox
                      << ", sourceFieldBox: " << sourceFieldBox << ", box: " << box
                      << ", intersectionBox: " << intersectionBox << std::endl;


            if constexpr (dimension == 1)
            {
                int iStartX = intersectionBox.lower(dirX);
                int iEndX   = intersectionBox.upper(dirX);

                for (int ix = iStartX; ix <= iEndX; ++ix)
                {
                    refiner(sourceField, destinationField, {{ix}});
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
                        refiner(sourceField, destinationField, {{ix, iy}});
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
                            refiner(sourceField, destinationField, {{ix, iy, iz}});
                        }
                    }
                }
            }
        }
    }
};
} // namespace PHARE::amr



#endif
