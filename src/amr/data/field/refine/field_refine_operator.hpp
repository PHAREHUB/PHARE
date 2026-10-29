#ifndef PHARE_FIELD_REFINE_OPERATOR_HPP
#define PHARE_FIELD_REFINE_OPERATOR_HPP

#include "core/def.hpp"
#include "core/data/vector.hpp"
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


// the kernel only ever sees plain fields: FieldT itself, or the tile field view of a tiled FieldT
template<typename FieldT>
struct refine_kernel_field
{
    using type = FieldT;
};
template<typename FieldT>
    requires(core::is_field_tile_set_v<FieldT>)
struct refine_kernel_field<FieldT>
{
    using type = FieldT::field_type;
};
template<typename FieldT>
using refine_kernel_field_t = refine_kernel_field<FieldT>::type;


// tiled refines borrow the untiled scratch field of the destination patch, whose id the
// messenger registers once the scratch is registered as a resource
struct RefineScratchUser
{
    virtual ~RefineScratchUser() = default;

    void registerScratchID(int const scratch_id) { scratch_id_ = scratch_id; }

protected:
    int scratch_id_ = -1;
};

void registerRefineScratchID(int const scratch_id, auto const&... ops)
{
    auto const set = [&](auto const& op) {
        if (auto const user = std::dynamic_pointer_cast<RefineScratchUser>(op))
            user->registerScratchID(scratch_id);
    };
    (set(ops), ...);
}


// tiled refines gather their coarse source into the untiled phare_scratch_field of the
// destination patch, or into fallback on the SAMRAI temporary patches (coarse interpolation
// levels) which have none
template<typename GridLayoutT, typename FieldT>
void with_refine_scratch(SAMRAI::hier::Patch const& patch, int const scratch_id, auto& fallback,
                         std::size_t const fallback_size, auto&& fn)
{
    if (scratch_id < 0)
        throw std::runtime_error("tiled field refine: registerScratchID was not called");
    if (patch.checkAllocated(scratch_id))
        fn(FieldData<GridLayoutT, typename FieldT::grid_type>::getField(patch, scratch_id));
    else
        fn(fallback.get_no_copy(fallback_size));
}


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

// Per overlap box: the coarse cells its stencils read are gathered into the scratch, the kernel
// fills each destination tile from it (the cells the tile owns, plus a ring so whole coarse cells
// round out as on the untiled patch), then the owners' values are copied into the other tiles'
// ghost copies, so inner tile ghosts agree
template<typename FieldT>
void refine_field(auto const& kernel, FieldT const& src, FieldT& dst, auto const& overlapBoxes,
                  auto const& centering, auto const& destFieldBox, auto const& sourceFieldBox,
                  auto const& ratio, auto& scratch)
    requires(core::is_field_tile_set_v<FieldT>)
{
    using Field_t              = FieldT::field_type;
    auto constexpr dim         = FieldT::dimension;
    auto const qty             = dst.physicalQuantity();
    auto const patch_field_box = phare_box_from<dim>(destFieldBox);
    auto const expander        = core::tile_owner_expander(dst);
    auto const& samrai_dim     = destFieldBox.getDim();
    auto const one             = SAMRAI::hier::IntVector::getOne(samrai_dim);
    auto const stencil         = SAMRAI::hier::IntVector{samrai_dim, kernel.coarseStencilWidth()};

    for (auto const& box : overlapBoxes)
    {
        auto const intersectionBox = destFieldBox * box;
        if (intersectionBox.empty())
            continue;

        auto const coarseBox = [&]() {
            auto coarse = intersectionBox;
            coarse.grow(one);
            coarse.coarsen(ratio);
            coarse.grow(stencil);
            return coarse * sourceFieldBox;
        }();
        auto const coarse_box = phare_box_from<dim>(coarseBox);
        if (coarse_box.size() > scratch.size())
            throw std::runtime_error("tiled field refine: scratch field too small");

        Field_t coarse{src.physicalQuantity(), scratch.data(), *coarse_box.shape().as_unsigned()};
        core::gather_owned(src, coarse_box, coarse, coarse_box.lower);

        auto const fine_box = phare_box_from<dim>(intersectionBox);
        for (auto& tile : dst())
        {
            auto const tile_gb = tile.layout().AMRGhostBoxFor(qty);
            if (auto const near = grow(expander(*tile), 1) * tile_gb)
                if (auto const region = fine_box * *near)
                    kernel.refineBox(coarse, tile(), samrai_box_from(*region), centering,
                                     samrai_box_from(tile_gb), coarseBox, ratio);
        }

        if (auto const written = grow(fine_box, 1) * patch_field_box)
            core::sync_owned(dst, *written);
    }
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
class KernelFieldRefineOperator : public SAMRAI::hier::RefineOperator, public RefineScratchUser
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

        if constexpr (core::is_field_tile_set_v<FieldT>)
            with_refine_scratch<GridLayoutT, FieldT>(
                destination, scratch_id_, tmp, sourceFieldBox.size(), [&](auto& scratch) {
                    refine_field(*kernel_, sourceField, destinationField, overlapBoxes,
                                 destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio,
                                 scratch);
                });
        else
            refine_field(*kernel_, sourceField, destinationField, overlapBoxes,
                         destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio);
    }

private:
    std::unique_ptr<Kernel_t> kernel_;

    static inline core::MinimizingVector<typename FieldT::type> tmp; // SAMRAI temporary patches
};


/**
 * @brief Tensor-field (rank 1 = vector) kernel operator. Per-component dispatch to the
 * same runtime kernel.
 */
template<typename TensorFieldData_t>
class KernelTensorFieldRefineOperator : public SAMRAI::hier::RefineOperator,
                                        public RefineScratchUser
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

            if constexpr (core::is_field_tile_set_v<GridT>)
                with_refine_scratch<GridLayoutT, GridT>(
                    destination, scratch_id_, tmp, sourceFieldBox.size(), [&](auto& scratch) {
                        refine_field(*kernel_, sourceFields[c], destinationFields[c], overlapBoxes,
                                     destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio,
                                     scratch);
                    });
            else
                refine_field(*kernel_, sourceFields[c], destinationFields[c], overlapBoxes,
                             destLayout.centering(qty), destFieldBox, sourceFieldBox, ratio);
        }
    }

private:
    std::unique_ptr<Kernel_t> kernel_;

    static inline core::MinimizingVector<typename GridT::type> tmp; // SAMRAI temporary patches
};


template<typename VectorFieldDataT>
using KernelVecFieldRefineOperator = KernelTensorFieldRefineOperator<VectorFieldDataT>;


} // namespace PHARE::amr

#endif
