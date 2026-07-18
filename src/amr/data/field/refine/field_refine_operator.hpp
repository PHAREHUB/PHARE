#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/field/field_box.hpp"

#include "amr/amr_constants.hpp"
#include "amr/data/field/field_data.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/resources_manager/tensor_field_resource.hpp"

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

    FieldRefineOperator()
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
            auto const fine_overlap = phare_box_from<dimension>(destFieldBox * box);
            core::FieldBox fine{destinationField, destLayout, fine_overlap};
            core::FieldBox coarse{sourceField, srcLayout, coarsen_box(fine_overlap)};
            core::operator_across_adjacent_levels(coarse, fine, refinementRatio, refiner);
        }
    }
};


template<std::size_t rank, typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
class TensorFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = GridLayoutT::implT;

    using Quantity         = extract_quantity_type<typename FieldT::physical_quantity_type>::type;
    using TensorFieldDataT = TensorFieldData<rank, GridLayoutT, FieldT, Quantity>;
    using TensorFieldOverlap_t = TensorFieldOverlap<rank>;

    static constexpr std::size_t N = TensorFieldDataT::N;

    TensorFieldRefineOperator()
        : SAMRAI::hier::RefineOperator{"TensorFieldRefineOperator"}

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
        auto const& destinationTensorFieldOverlap
            = dynamic_cast<TensorFieldOverlap_t const&>(destinationOverlap);
        auto const& srcData      = source.getPatchData(sourceId);
        auto const& destData     = destination.getPatchData(destinationId);
        auto& destinationFields  = TensorFieldDataT::getFields(destination, destinationId);
        auto const& destLayout   = TensorFieldDataT::getLayout(destination, destinationId);
        auto const& sourceFields = TensorFieldDataT::getFields(source, sourceId);
        auto const& srcLayout    = TensorFieldDataT::getLayout(source, sourceId);

        PHARE_LOG_SCOPE(2, "TensorFieldRefineOperator::refine::" + sourceFields[0].name());

        // We assume that quantity are all the same.
        // Note that an assertion will be raised in refineIt operator
        for (std::uint16_t c = 0; c < N; ++c)
        {
            auto const& overlapBoxes
                = destinationTensorFieldOverlap[c]->getDestinationBoxContainer();
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
                auto const fine_overlap = phare_box_from<dimension>(destFieldBox * box);
                core::FieldBox fine{destinationFields[c], destLayout, fine_overlap};

                auto const coarse_overlap
                    = *(coarsen_box(fine_overlap) * srcLayout.AMRGhostBoxFor(sourceFields[c]));
                core::FieldBox coarse{sourceFields[c], srcLayout, coarse_overlap};
                core::operator_across_adjacent_levels(coarse, fine, refinementRatio, refiner);
            }
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
using VecFieldRefineOperator
    = TensorFieldRefineOperator</*rank=*/1, GridLayoutT, FieldT, FieldRefinerPolicy>;


} // namespace PHARE::amr



#endif
