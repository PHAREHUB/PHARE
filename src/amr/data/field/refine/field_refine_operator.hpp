#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "amr/data/field/field_data.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/resources_manager/tensor_field_resource.hpp"

// #include "field_linear_refine.hpp"
#include "field_refiner.hpp"

#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/RefineOperator.h>


#include <cstddef>
#include <stdexcept>


namespace PHARE::amr
{


template<typename Dst>
void refine_field(Dst& destinationField, auto& sourceField, auto& intersectionBox, auto& refiner)
{
    for (auto const bix : phare_box_from<Dst::dimension>(intersectionBox))
        refiner(sourceField, destinationField, bix);
}

template<typename Refiner, typename FieldT>
void refine_field(FieldT& dst, auto const& dstBox, auto const& dstLayout, FieldT const& src,
                  auto const& srcBox, auto& overlap, auto ratio)
// requires(not core::is_field_tile_set_v<FieldT>)
{
    auto const& qty = dst.physicalQuantity();

    Refiner refiner{dstLayout.centering(qty), dstBox, srcBox, ratio};

    for (auto const& box : overlap.getDestinationBoxContainer())
    {
        // we compute the intersection with the destination,
        // and then we apply the refine operation on each fine index.
        auto intersectionBox = dstBox * box;
        refine_field(dst, src, intersectionBox, refiner);
    }
}

// template<typename Refiner, typename FieldT>
// void refine_field(FieldT& dst, auto const& dstBox, auto const& dstLayout, FieldT const& src,
//                   auto const& srcBox, auto& overlap, auto ratio)
//     requires(core::is_field_tile_set_v<FieldT>)
// {
//     auto const& qty = dst.physicalQuantity();

//     throw std::runtime_error("finish");
// }



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

        refine_field<FieldRefinerPolicy>(destinationField, destFieldBox, destLayout, sourceField,
                                         sourceFieldBox, destinationFieldOverlap, ratio);
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

            refine_field<FieldRefinerPolicy>(destinationFields[c], destFieldBox, destLayout,
                                             sourceFields[c], sourceFieldBox,
                                             *destinationTensorFieldOverlap[c], ratio);
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename FieldRefinerPolicy>
using VecFieldRefineOperator
    = TensorFieldRefineOperator</*rank=*/1, GridLayoutT, FieldT, FieldRefinerPolicy>;


} // namespace PHARE::amr



#endif
