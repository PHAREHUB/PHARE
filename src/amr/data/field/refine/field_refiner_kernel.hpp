#ifndef PHARE_FIELD_REFINER_KERNEL_HPP
#define PHARE_FIELD_REFINER_KERNEL_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>

#include <array>
#include <cstddef>
#include <memory>


namespace PHARE::amr
{

/**
 * @brief Supported field-refinement orders, shared by RefinementConfig and the kernel factories.
 *
 * Linear (order 2) is the only supported value for now; validating a raw dict int against this
 * enum happens in exactly one place (RefinementConfig::FROM), so the factories below no longer
 * need their own runtime order check.
 */
enum class FieldRefinementOrder { Linear = 2 };

/**
 * @brief Runtime-dispatched field-refinement seam.
 *
 * A kernel is constructed once per refine operator and applied per overlap box: refineBox()
 * receives one intersection box, together with {centering, destFieldBox, sourceFieldBox, ratio},
 * and loops over its fine indices internally, so the virtual-dispatch cost is paid once per box
 * rather than once per index.
 *
 * Concrete kernels (composite Linear, magnetic shared-face) implement refineBox().
 */
template<typename GridLayoutT, typename FieldT>
struct IFieldRefineKernel
{
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    virtual void refineBox(FieldT const& sourceField, FieldT& destinationField,
                           SAMRAI::hier::Box const& intersectionBox,
                           std::array<core::QtyCentering, dimension> const& centering,
                           SAMRAI::hier::Box const& destFieldBox,
                           SAMRAI::hier::Box const& sourceFieldBox,
                           SAMRAI::hier::IntVector const& ratio) const
        = 0;

    /**
     * @brief Coarse-cell stencil half-width this kernel reads around each anchor.
     *
     * SAMRAI provisions coarse (source) ghost layers from RefineOperator::getStencilWidth before
     * prolongation. order 2 reads ±1 coarse cell (both the dual ±¼ ladder and the primal
     * midpoint), so it is order/2. Reported up through the holding operator.
     */
    virtual int coarseStencilWidth() const = 0;

    virtual ~IFieldRefineKernel() = default;
};


/**
 * @brief Build a composite field-refinement kernel for a given order.
 *
 * order: 2 = Linear (per the dual ±1/4 ladder). Definition lives with the concrete composite
 * kernels (composite_field_refiner.hpp); declared here so the additive operators and the
 * messengers can depend only on the seam.
 */
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeRefineKernel(FieldRefinementOrder order);

/**
 * @brief Build the stage-1 magnetic refinement kernel of the ADPT div-free prolongation.
 *
 * Fills ALL fine faces per component (shared and interior) with the composite tensor stencils;
 * the stage-2 partner ADPTMagneticRefinePatchStrategy::postprocessRefine then applies the
 * order-independent divB touch-up. Stage 1 by itself makes no ∇·B claim — it only reproduces
 * each component's coarse values on the fine grid via the tensor stencils above; ∇·B-freeness of
 * the composite result is established solely by the stage-2 touch-up, which equalizes the 2^d
 * fine subzone divergences of each coarse zone once stage 1 has filled every face (see
 * ADPTMagneticRefinePatchStrategy's class docs).
 */
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeMagneticRefineKernel(FieldRefinementOrder order);


} // namespace PHARE::amr


#endif
