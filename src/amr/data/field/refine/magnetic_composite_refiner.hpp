#ifndef PHARE_MAGNETIC_COMPOSITE_REFINER_HPP
#define PHARE_MAGNETIC_COMPOSITE_REFINER_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"

#include "coarse_cell_round_out.hpp"
#include "composite_field_refiner.hpp"
#include "field_refiner_kernel.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>


namespace PHARE::amr
{

/**
 * @brief Stage 1 of the Balsara ADPT divB-free B prolongation: fill every fine face of a B
 * component from its own coarse faces.
 *
 * Per component the fill is the composite one (composite_field_refiner.hpp): the primal-even
 * direction (collocated with coarse) an exact copy, primal-odd (a new fine face) a half-point
 * interpolation, the dual directions the ±¼ child ladder. This stage makes no divB claim — stage
 * 2, the cross-component touch-up (adpt_magnetic_refine_patch_strategy.hpp), is what establishes
 * it.
 *
 * What is magnetic here, and the whole reason this is a kernel of its own rather than the plain
 * composite one, is the region gathered: stage 2 reads the shared faces of a whole coarse cell, so
 * stage 1 has to have written them.
 */
template<typename GridLayoutT, typename FieldT, std::size_t order>
class MagneticCompositeRefiner final : public IFieldRefineKernel<GridLayoutT, FieldT>
{
    static constexpr std::size_t dimension = GridLayoutT::dimension;

    using Fill = detail::CompositeFill<GridLayoutT, FieldT, order>;

public:
    void refineBox(FieldT const& sourceField, FieldT& destinationField,
                   SAMRAI::hier::Box const& intersectionBox,
                   std::array<core::QtyCentering, dimension> const& centering,
                   SAMRAI::hier::Box const& destFieldBox, SAMRAI::hier::Box const& sourceFieldBox,
                   SAMRAI::hier::IntVector const& ratio) const override
    {
        Fill::fill(sourceField, destinationField,
                   gatherRegion_(intersectionBox, centering, destFieldBox), centering, destFieldBox,
                   sourceFieldBox, ratio);
    }

    int coarseStencilWidth() const override { return Fill::coarseStencilWidth(); }

private:
    /**
     * @brief The fine region to gather for an overlap: whole coarse cells, so the stage-2 touch-up
     * finds every face it reads already written.
     *
     * A B component is primal in at most one direction, its face normal. Without one — the
     * out-of-plane components, By and Bz in 1D, Bz in 2D — the component has no interior face,
     * hence no divB neighbourhood to protect and nothing to round out.
     *
     * With one, round out to whole coarse cells (coarse_cell_round_out.hpp), then re-clip: the
     * overlap arrives already clipped to the allocation, so rounding can push it back out. Widening
     * the gather is safe because the fill writes only still-NaN fine indices, and the fine indices
     * rounded in share the coarse anchors, and the stencil reach, of those already gathered.
     */
    static SAMRAI::hier::Box
    gatherRegion_(SAMRAI::hier::Box const& overlap,
                  std::array<core::QtyCentering, dimension> const& centering,
                  SAMRAI::hier::Box const& destFieldBox)
    {
        auto const primalCount
            = std::count_if(centering.begin(), centering.end(),
                            [](auto const c) { return c == core::QtyCentering::primal; });

        if (primalCount > 1)
            throw std::runtime_error(
                "magnetic refiner expects at most one primal (normal) direction");

        if (primalCount == 0)
            return overlap;

        return roundFieldBoxOutToCoarseCells<dimension>(overlap, centering) * destFieldBox;
    }
};


// Single-enumerator dispatch: the order is validated once in RefinementConfig::FROM.
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeMagneticRefineKernel([[maybe_unused]] FieldRefinementOrder const order)
{
    return std::make_unique<MagneticCompositeRefiner<
        GridLayoutT, FieldT, static_cast<std::size_t>(FieldRefinementOrder::Linear)>>();
}


} // namespace PHARE::amr


#endif
