#ifndef PHARE_FIELD_REFINER_KERNEL_HPP
#define PHARE_FIELD_REFINER_KERNEL_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

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
 * Linear (order 2) is the only supported value for now. A raw dict int is validated against this
 * enum in exactly one place, RefinementConfig::FROM; the factories below take the enum and so need
 * no order check of their own.
 */
enum class FieldRefinementOrder { Linear = 2 };

/**
 * @brief Runtime-dispatched field-refinement seam.
 *
 * The refinement order is a runtime dict value, so it must not become a template parameter of the
 * refine operators and messengers: that would carry it into the build permutation id and multiply
 * the compiled permutations. This seam is what keeps it runtime — operators and messengers depend
 * on the interface, never on the stencil tables.
 *
 * refineBox() takes one overlap box and loops over its fine indices internally, so the virtual
 * call is paid once per box rather than once per index.
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
 * Defined with the concrete kernels (composite_field_refiner.hpp); declared here so the additive
 * operators and the messengers depend only on the seam.
 */
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeRefineKernel(FieldRefinementOrder order);

/**
 * @brief Build the stage-1 magnetic refinement kernel of the ADPT div-free prolongation.
 *
 * Fills every fine face per component with the composite tensor stencils, and makes no ∇·B claim:
 * ∇·B-freeness comes from the stage-2 touch-up (adpt_magnetic_refine_patch_strategy.hpp).
 */
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeMagneticRefineKernel(FieldRefinementOrder order);


} // namespace PHARE::amr


#endif
