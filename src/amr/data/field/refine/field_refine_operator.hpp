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



/**
 * @brief Runtime-dispatched field refine operator.
 *
 * Holds a unique_ptr<IFieldRefineKernel> chosen at construction (makeRefineKernel(order)) and
 * forwards each overlap box to the kernel, which decides the stencil from the refinement order.
 * Nothing is ever shared: the factory mints a fresh kernel per operator, and the kernels are
 * stateless (their stencil tables are static, held behind the kernel's vtable).
 */
template<typename GridLayoutT, typename FieldT>
class KernelFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using PhysicalQuantity                 = FieldT::physical_quantity_type;
    using FieldDataT                       = FieldData<GridLayoutT, FieldT>;
    using Kernel_t                         = IFieldRefineKernel<GridLayoutT, FieldT>;

    explicit KernelFieldRefineOperator(std::unique_ptr<Kernel_t> kernel)
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
        auto const& overlapBoxes            = destinationFieldOverlap.getDestinationBoxContainer();

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
    std::unique_ptr<Kernel_t> kernel_;
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

    explicit KernelTensorFieldRefineOperator(std::unique_ptr<Kernel_t> kernel)
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
    std::unique_ptr<Kernel_t> kernel_;
};


template<typename VectorFieldDataT>
using KernelVecFieldRefineOperator = KernelTensorFieldRefineOperator<VectorFieldDataT>;


} // namespace PHARE::amr



#endif
