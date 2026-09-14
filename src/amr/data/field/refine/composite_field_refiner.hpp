#ifndef PHARE_COMPOSITE_FIELD_REFINER_HPP
#define PHARE_COMPOSITE_FIELD_REFINER_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "coarse_cell_round_out.hpp"
#include "field_refiner_kernel.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>


namespace PHARE::amr
{

/**
 * @brief Runtime field-refinement kernel built from the two 1-D primitives tensor-producted.
 *
 * For ratio 2, each fine index maps to a coarse anchor I = floor(f/2) (toCoarseIndex) and a parity
 * p = f − 2I ∈ {0,1} per direction. Per direction the 1-D weight row is chosen by centering+parity:
 *   - primal, p=0 : exact copy (coincident node)            → {offset 0, 1.0}
 *   - primal, p=1 : half-point midpoint (Primitive A.2)     → directionalInterp<dir,PrimalToDual>
 *   - dual,   p   : ±¼ child ladder (Primitive B)           → directionalProlongation<dir, σ=2p−1>
 * (both at <order>). The multi-D stencil is the outer product of the rows (consteval
 * tensorProduct), which sizes each multi-D stencil exactly from its rows' own lengths (copy=1,
 * half-point=2, dual-σ±=3) — no padding to a common length; the full set of centering×parity
 * combinations (4^dim of them) is a compile-time table, indexed at runtime per fine index by
 * centering+parity. Offsets are relative to the anchor I; values are gathered from the coarse
 * field.
 */
template<typename GridLayoutT, typename FieldT, std::size_t order, bool isMagnetic = false>
class CompositeFieldRefiner : public IFieldRefineKernel<GridLayoutT, FieldT>
{
    static constexpr std::size_t dimension = GridLayoutT::dimension;
    using GridLayoutImpl                   = typename GridLayoutT::implT;
    using Point_t                          = core::Point<int, dimension>;
    using WeightPoint_t                    = core::WeightPoint<dimension>;

    static_assert(order == 2, "composite refiner ladder is order 2 (Linear)");

public:
    void refineBox(FieldT const& sourceField, FieldT& destinationField,
                   SAMRAI::hier::Box const& intersectionBox,
                   std::array<core::QtyCentering, dimension> const& centering,
                   SAMRAI::hier::Box const& destFieldBox, SAMRAI::hier::Box const& sourceFieldBox,
                   SAMRAI::hier::IntVector const& ratio) const override
    {
        for (std::size_t d = 0; d < dimension; ++d)
            if (ratio(d) != 2)
                throw std::runtime_error("CompositeFieldRefiner supports refinement ratio 2 only");

        // Stage 1 fills every fine face, interior ones included; the stage-2 touch-up then
        // overwrites the interior faces. What matters here is only whether the component has a
        // normal (face) direction at all, since that gates the round-out below. A B component is
        // primal in at most one direction; an out-of-plane component (By,Bz in 1D, Bz in 2D) has
        // none, hence no interior face and no ∇·B neighbourhood to protect.
        [[maybe_unused]] bool hasNormal = false;
        if constexpr (isMagnetic)
        {
            auto const primalCount
                = std::count_if(centering.begin(), centering.end(),
                                [](auto const c) { return c == core::QtyCentering::primal; });
            if (primalCount > 1)
                throw std::runtime_error(
                    "magnetic refiner expects at most one primal (normal) direction");
            hasNormal = (primalCount == 1);
        }

        // The stage-2 touch-up reads the shared faces of a whole coarse cell, so this gather must
        // satisfy the same whole-coarse-cell invariant (coarse_cell_round_out.hpp). Round out,
        // then re-clip: the box arrives already clipped, so rounding can push it past the
        // allocation. Safe — assignFine_ writes only still-NaN slots, and the rounded-in fine
        // indices share the coarse anchors, and stencil reach, of those already gathered.
        auto const region = [&]() -> SAMRAI::hier::Box {
            if constexpr (isMagnetic)
                if (hasNormal)
                    return roundFieldBoxOutToCoarseCells<dimension>(intersectionBox, centering)
                           * destFieldBox;
            return intersectionBox;
        }();

        // Centering is fixed over the box, parity varies per fine index. The 2-bits-per-direction
        // index built below (bit0 = parity, bit1 = dual) selects the matching compile-time stencil.
        for (auto const fineIndex : phare_box_from<dimension>(region))
        {
            auto const anchor = toCoarseIndex<dimension>(fineIndex);

            std::size_t combined = 0; // 2 bits/dir: (dual<<1)|parity
            for (std::size_t d = 0; d < dimension; ++d)
            {
                auto const parity         = static_cast<std::size_t>(fineIndex[d] - 2 * anchor[d]);
                std::size_t const dualBit = (centering[d] == core::QtyCentering::dual) ? 1u : 0u;
                combined |= ((dualBit << 1) | parity) << (2u * d);
            }

            auto const anchorLocal = AMRToLocal(anchor, sourceFieldBox);
            auto const fineLocal   = AMRToLocal(fineIndex, destFieldBox);

            double const value = gatherers_[combined](sourceField, anchorLocal);

            assignFine_(destinationField, fineLocal, value);
        }
    }

    // order 2 reads ±1 coarse cell (max |offset| over both 1-D primitives).
    int coarseStencilWidth() const override { return order / 2; }

private:
    // 2 bits of choice per direction (dual<<1 | parity) → 4^dim combinations.
    static constexpr std::size_t nCombined = std::size_t{1} << (2 * dimension);

    // 1-D weight row for direction d and the compile-time choice (dual<<1 | parity).
    template<std::size_t d, std::size_t choice>
    static consteval auto oneDRow_()
    {
        if constexpr (choice == 0) // primal, parity 0: coincident node, exact copy
            return std::array{WeightPoint_t{Point_t{}, 1.0}};
        else if constexpr (choice == 1) // primal, parity 1: half-point midpoint
            return GridLayoutImpl::template directionalInterp<
                d, GridLayoutImpl::InterpDir::PrimalToDual, order>();
        else if constexpr (choice == 2) // dual, parity 0: σ = −1 child
            return GridLayoutImpl::template directionalProlongation<d, -1, order>();
        else // choice == 3 : dual, parity 1: σ = +1 child
            return GridLayoutImpl::template directionalProlongation<d, +1, order>();
    }

    // full multi-D stencil (exact size) for the combined choice index (2 bits per direction).
    template<std::size_t combined>
    static consteval auto makeStencil_()
    {
        if constexpr (dimension == 1)
            return oneDRow_<0, combined & 3u>();
        else if constexpr (dimension == 2)
            return GridLayoutImpl::template tensorProduct<0, 1>(
                oneDRow_<0, combined & 3u>(), oneDRow_<1, (combined >> 2) & 3u>());
        else
            return GridLayoutImpl::template tensorProduct<0, 1, 2>(
                oneDRow_<0, combined & 3u>(), oneDRow_<1, (combined >> 2) & 3u>(),
                oneDRow_<2, (combined >> 4) & 3u>());
    }

    // gather the coarse contributions for one fine index using the exact stencil of choice
    // `combined` (its size is baked in at compile time — no zero-coef points iterated).
    template<std::size_t combined>
    static double gatherStencil_(FieldT const& src, Point_t const& anchorLocal)
    {
        double value = 0.;
        for (auto const& w : makeStencil_<combined>())
            value += w.coef * src(anchorLocal + w.indexes);
        return value;
    }

    // runtime centering+parity → compile-time stencil: one gather fn per combination, dispatched
    // by a plain index.
    using Gatherer_t                 = double (*)(FieldT const&, Point_t const&);
    static constexpr auto gatherers_ = core::for_N_make_array<nCombined>(
        [](auto ic) { return Gatherer_t{&gatherStencil_<ic()>}; });

    // preserve the legacy NaN-guard: only fill fine indices not already set
    static void assignFine_(FieldT& dst, Point_t const& fineLocal, double const value)
    {
        if (auto& dst_val = dst(fineLocal); std::isnan(dst_val))
            dst_val = value;
    }
};


// ---- factory (declared in field_refiner_kernel.hpp) ---------------------------------------------

// FieldRefinementOrder has a single enumerator and RefinementConfig::FROM is the only place a
// dict value is validated into it, so there is nothing to branch on. A second order adds a switch
// here, one case per compile-time stencil.
template<typename GridLayoutT, typename FieldT>
std::unique_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeRefineKernel([[maybe_unused]] FieldRefinementOrder const order)
{
    return std::make_unique<CompositeFieldRefiner<
        GridLayoutT, FieldT, static_cast<std::size_t>(FieldRefinementOrder::Linear)>>();
}


} // namespace PHARE::amr


#endif
