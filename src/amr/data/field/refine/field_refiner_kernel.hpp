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
 * order: 2 = Linear (per the dual ±1/4 ladder). Definition lives with the concrete
 * composite kernels (Step 3); declared here so the additive operators and the messengers can
 * depend only on the seam.
 */
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>> makeRefineKernel(int order);

/**
 * @brief Build the stage-1 magnetic refinement kernel of the ADPT div-free prolongation.
 *
 * Fills ALL fine faces per component (shared and interior) with the composite tensor stencils;
 * the stage-2 partner ADPTMagneticRefinePatchStrategy::postprocessRefine then applies the
 * order-independent divB touch-up. The shared-face tangential correction is antisymmetric in the
 * child sign, so ∇·B is preserved before the touch-up ever runs.
 */
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>> makeMagneticRefineKernel(int order);


} // namespace PHARE::amr


#endif
