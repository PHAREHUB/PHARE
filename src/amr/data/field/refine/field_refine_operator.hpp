#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP

#include "core/def.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/field/field_box.hpp"

#include "amr/data/field/field_data.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"

#include "field_refiner_kernel.hpp"

#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/hier/RefineOperator.h>

#include <cstddef>
#include <type_traits>

namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;


// the kernel only ever sees plain fields: FieldT itself, or the grid a tiled FieldT reduces to
template<typename FieldT>
using refine_kernel_field_t
    = std::decay_t<decltype(core::reduce_single(std::declval<FieldT const&>()))>;


template<typename FieldT>
void refine_field(auto const& kernel, FieldT const& src, FieldT& dst, auto const& overlapBoxes,
                  auto const& centering, auto const& destFieldBox, auto const& sourceFieldBox,
                  auto const& ratio)
    requires(not core::is_field_tile_set_v<FieldT>)
{
    for (auto const& box : overlapBoxes)
    {
        auto const intersectionBox = destFieldBox * box;
        kernel.refineBox(src, dst, intersectionBox, centering, destFieldBox, sourceFieldBox, ratio);
    }
}

// TODO: no dynamic allocation -- borrow a pre-declared untiled temporary instead of reducing
// into fresh grids (SAMRAI's coarse-interpolation temporary patches have none allocated)
template<typename FieldT>
void refine_field(auto const& kernel, FieldT const& src, FieldT& dst, auto const& overlapBoxes,
                  auto const& centering, auto const& destFieldBox, auto const& sourceFieldBox,
                  auto const& ratio)
    requires(core::is_field_tile_set_v<FieldT>)
{
    auto const flat_src = core::reduce_single(src);
    auto flat_dst       = core::reduce_single(dst); // keeps NaNs, the kernel only fills those
    refine_field(kernel, flat_src, flat_dst, overlapBoxes, centering, destFieldBox, sourceFieldBox,
                 ratio);
    core::copy_fields(dst, *flat_dst); // every tile's ghost box, so inner tile ghosts agree
}



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
    using Kernel_t = IFieldRefineKernel<GridLayoutT, refine_kernel_field_t<FieldT>>;

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

        refine_field(*kernel_, sourceField, destinationField, overlapBoxes,
                     destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio);
    }

private:
    std::unique_ptr<Kernel_t> kernel_;
};


/**
 * @brief Tensor-field (rank 1 = vector) kernel operator. Per-component dispatch to the
 * same runtime kernel.
 */
template<typename TensorFieldData_t>
class KernelTensorFieldRefineOperator : public SAMRAI::hier::RefineOperator
{
public:
    using GridLayoutT                      = TensorFieldData_t::gridlayout_type;
    using GridT                            = TensorFieldData_t::grid_type;
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using TensorFieldDataT     = TensorFieldData_t;
    using TensorFieldOverlap_t = TensorFieldOverlap<TensorFieldData_t::rank>;
    using Kernel_t             = IFieldRefineKernel<GridLayoutT, refine_kernel_field_t<GridT>>;

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

            refine_field(*kernel_, sourceFields[c], destinationFields[c], overlapBoxes,
                         destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio);
        }
    }

private:
    std::unique_ptr<Kernel_t> kernel_;
};


template<typename VectorFieldDataT>
using KernelVecFieldRefineOperator = KernelTensorFieldRefineOperator<VectorFieldDataT>;


} // namespace PHARE::amr

#endif
