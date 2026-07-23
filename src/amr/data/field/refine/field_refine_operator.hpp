#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "amr/data/field/field_data.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"

#include "field_refiner_kernel.hpp"

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
    for (auto const bix : phare_box_from<Dst::dimension>(intersectionBox))
        refiner(sourceField, destinationField, bix);
}


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
            auto intersectionBox = destFieldBox * box;
            refine_field(destinationField, sourceField, intersectionBox, refiner);
        }
    }
};


template<typename TensorFieldData_t, typename FieldRefinerPolicy>
class TensorFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    using GridLayoutT                      = TensorFieldData_t::gridlayout_type;
    using FieldT                           = TensorFieldData_t::grid_type::field_type;
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = GridLayoutT::implT;

    using TensorFieldDataT     = TensorFieldData_t;
    using TensorFieldOverlap_t = TensorFieldOverlap<TensorFieldData_t::rank>;

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

            FieldRefinerPolicy refiner{destLayout.centering(qty), destFieldBox, sourceFieldBox,
                                       ratio};

            for (auto const& box : overlapBoxes)
            {
                // we compute the intersection with the destination,
                // and then we apply the refine operation on each fine index.
                auto const intersectionBox = destFieldBox * box;
                refine_field(destinationFields[c], sourceFields[c], intersectionBox, refiner);
            }
        }
    }
};

template<typename VectorFieldDataT, typename FieldRefinerPolicy>
using VecFieldRefineOperator = TensorFieldRefineOperator<VectorFieldDataT, FieldRefinerPolicy>;

/**
 * @brief Additive, runtime-dispatched counterpart of FieldRefineOperator.
 *
 * Holds a shared_ptr<IFieldRefineKernel> chosen at construction (makeRefineKernel(order,limiter))
 * and forwards each overlap box to the kernel. The legacy templated FieldRefineOperator above is
 * left byte-identical; this operator is used ONLY on the composite order=2/4 path, so the default
 * (order absent) carries no behavior change.
 */
template<typename GridLayoutT, typename FieldT>
class KernelFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using PhysicalQuantity                 = typename FieldT::physical_quantity_type;
    using FieldDataT                       = FieldData<GridLayoutT, FieldT>;
    using Kernel_t                         = IFieldRefineKernel<GridLayoutT, FieldT>;

    explicit KernelFieldRefineOperator(std::shared_ptr<Kernel_t> kernel)
        : SAMRAI::hier::RefineOperator{"KernelFieldRefineOperator"}
        , kernel_{std::move(kernel)}
    {
    }

    virtual ~KernelFieldRefineOperator() = default;

    NO_DISCARD int getOperatorPriority() const override { return 0; }

    NO_DISCARD SAMRAI::hier::IntVector
    getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector(dim, kernel_->coarseStencilWidth());
    }

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

        auto const& qty     = destinationField.physicalQuantity();
        auto const destData = destination.getPatchData(destinationId);
        auto const srcData  = source.getPatchData(sourceId);

        auto const destFieldBox
            = FieldGeometry::toFieldBox(destData->getGhostBox(), qty, destLayout);
        auto const sourceFieldBox
            = FieldGeometry::toFieldBox(srcData->getGhostBox(), qty, srcLayout);

        auto const centering = destLayout.centering(qty);

        for (auto const& box : overlapBoxes)
        {
            auto const intersectionBox = destFieldBox * box;
            kernel_->refineBox(sourceField, destinationField, intersectionBox, centering,
                               destFieldBox, sourceFieldBox, ratio);
        }
    }

private:
    std::shared_ptr<Kernel_t> kernel_;
};


/**
 * @brief Tensor-field (rank 1 = vector) additive kernel operator. Per-component dispatch to the
 * same runtime kernel.
 */
template<typename TensorFieldData_t>
class KernelTensorFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    using GridLayoutT                      = TensorFieldData_t::gridlayout_type;
    using FieldT                           = TensorFieldData_t::grid_type;
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using TensorFieldDataT     = TensorFieldData_t;
    using TensorFieldOverlap_t = TensorFieldOverlap<TensorFieldData_t::rank>;
    using Kernel_t             = IFieldRefineKernel<GridLayoutT, FieldT>;

    static constexpr std::size_t N = TensorFieldDataT::N;

    explicit KernelTensorFieldRefineOperator(std::shared_ptr<Kernel_t> kernel)
        : SAMRAI::hier::RefineOperator{"KernelTensorFieldRefineOperator"}
        , kernel_{std::move(kernel)}
    {
    }

    virtual ~KernelTensorFieldRefineOperator() = default;

    NO_DISCARD int getOperatorPriority() const override { return 0; }

    NO_DISCARD SAMRAI::hier::IntVector
    getStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector(dim, kernel_->coarseStencilWidth());
    }

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

            auto const centering = destLayout.centering(qty);

            for (auto const& box : overlapBoxes)
            {
                auto const intersectionBox = destFieldBox * box;
                kernel_->refineBox(sourceFields[c], destinationFields[c], intersectionBox,
                                   centering, destFieldBox, sourceFieldBox, ratio);
            }
        }
    }

private:
    std::shared_ptr<Kernel_t> kernel_;
};


template<typename VectorFieldDataT>
using KernelVecFieldRefineOperator = KernelTensorFieldRefineOperator<VectorFieldDataT>;


} // namespace PHARE::amr



#endif
