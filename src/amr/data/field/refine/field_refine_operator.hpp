#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP


#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "core/def/phare_mpi.hpp"

#include "core/def.hpp"
#include "amr/data/field/field_data.hpp"

#include "core/hybrid/hybrid_quantities.hpp"
#include "field_linear_refine.hpp"

#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/RefineOperator.h>


#include <cstddef>


namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;



template<typename Dst>
void refine_field(Dst& destinationField, auto& sourceField, auto& intersectionBox, auto& refiner)
{
    auto constexpr static dimension = Dst::dimension;

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
        // Note that an assertion will be raised in refineIt operator
        auto const& qty     = destinationField.physicalQuantity();
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
            // and then we apply the refine operation on each fine index.
            auto intersectionBox = destFieldBox * box;
            refine_field(destinationField, sourceField, intersectionBox, refiner);
        }
    }
};


            std::cout << "debug refine operator: "
                      << "destinationFieldBox: " << destFieldBox
                      << ", sourceFieldBox: " << sourceFieldBox << ", box: " << box
                      << ", intersectionBox: " << intersectionBox << std::endl;


template<std::size_t rank, typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
class TensorFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = GridLayoutT::implT;

    using TensorFieldDataT = TensorFieldData<rank, GridLayoutT, FieldT, core::HybridQuantity>;

    static constexpr std::size_t N = TensorFieldDataT::N;

    TensorFieldRefineOperator(bool node_only = false)
        : SAMRAI::hier::RefineOperator{"FieldRefineOperator"}

    {
    }

    virtual ~TensorFieldRefineOperator() = default;

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
        auto const& destinationFieldOverlap = dynamic_cast<FieldOverlap const&>(destinationOverlap);
        auto const& overlapBoxes            = destinationFieldOverlap.getDestinationBoxContainer();
        auto const& srcData                 = source.getPatchData(sourceId);
        auto const& destData                = destination.getPatchData(destinationId);
        auto& destinationFields  = TensorFieldDataT::getFields(destination, destinationId);
        auto const& destLayout   = TensorFieldDataT::getLayout(destination, destinationId);
        auto const& sourceFields = TensorFieldDataT::getFields(source, sourceId);
        auto const& srcLayout    = TensorFieldDataT::getLayout(source, sourceId);

        // We assume that quantity are all the same.
        // Note that an assertion will be raised in refineIt operator
        for (std::uint16_t c = 0; c < N; ++c)
        {
            auto const& qty     = destinationFields[c].physicalQuantity();
            using FieldGeometry = FieldGeometry<GridLayoutT, std::decay_t<decltype(qty)>>;

            auto const destFieldBox
                = FieldGeometry::toFieldBox(destData->getGhostBox(), qty, destLayout);
            auto const sourceFieldBox
                = FieldGeometry::toFieldBox(srcData->getGhostBox(), qty, srcLayout);

            FieldRefinerPolicy refiner{destLayout.centering(qty), destFieldBox, sourceFieldBox,
                                       ratio};

            for (auto const& box : overlapBoxes)
            {
                // we compute the intersection with the destination,
                // and then we apply the refine operation on each fine index.
                auto intersectionBox = destFieldBox * box;
                refine_field(destinationFields[c], sourceFields[c], intersectionBox, refiner);
            }
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
using VecFieldRefineOperator
    = TensorFieldRefineOperator<1, GridLayoutT, FieldT, FieldRefinerPolicy>;


} // namespace PHARE::amr



#endif
