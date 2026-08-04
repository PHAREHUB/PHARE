#ifndef PHARE_COMPOSITE_FIELD_REFINER_HPP
#define PHARE_COMPOSITE_FIELD_REFINER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "coarse_cell_round_out.hpp"
#include "field_refiner_kernel.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>

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
 * tensorProduct). Each 1-D row is padded to length 3 with zero-coef points so all stencils share
 * one type; the full set of centering×parity combinations (4^dim of them) is a compile-time
 * table, indexed at runtime per fine index by centering+parity. Offsets are relative to the anchor
 * I; values are gathered from the coarse field.
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

        // Stage 1 of the Balsara ADPT div-free prolongation fills EVERY fine face of a B
        // component (shared and interior alike) from its coarse faces; ownership of divB-safety
        // for the interior (odd-normal) faces moves to the stage-2 touch-up
        // (ADPTMagneticRefinePatchStrategy::postprocessRefine), which reads them back and
        // overwrites them with a divergence-equalizing correction. We still need to know whether
        // this component HAS a normal (face) direction, though: that alone gates the whole-
        // coarse-cell round-out below, which the stage-2 touch-up depends on (see class notes for
        // the "hasNormal" split). A B component is primal in at most one direction (its face
        // normal); an out-of-plane / reduced-dimension component (By,Bz in 1D, Bz in 2D) has none,
        // has no interior face and no ∇·B neighbourhood to protect, so it needs no rounding.
        [[maybe_unused]] bool hasNormal = false;
        if constexpr (isMagnetic)
        {
            int primalCount = 0;
            for (std::size_t d = 0; d < dimension; ++d)
                if (centering[d] == core::QtyCentering::primal)
                    ++primalCount;
            if (primalCount > 1)
                throw std::runtime_error(
                    "magnetic refiner expects at most one primal (normal) direction");
            hasNormal = (primalCount == 1);
        }

        // The stage-2 touch-up that owns the interior faces reads the shared faces of the whole
        // coarse cell, so the shared faces of every coarse cell the region touches must exist —
        // the same whole-coarse-cell invariant the touch-up rounds its own region out to (see
        // coarse_cell_round_out.hpp and MagneticPatchStrategyBase::reconstructionRegion). Round
        // out here too, then re-clip: the box arrives already clipped to the destination, so
        // rounding can push it past the allocation. This cannot corrupt anything — assignFine_
        // writes only slots that are still NaN, so the extra layer of shared faces fills holes and
        // never overwrites valid data. It also needs no extra coarse data: for both centerings the
        // rounded-in fine indices share the coarse anchors, and the stencil reach around them, of
        // the indices already gathered.
        auto const region = [&]() -> SAMRAI::hier::Box {
            if constexpr (isMagnetic)
                if (hasNormal)
                    return roundFieldBoxOutToCoarseCells<dimension>(intersectionBox, centering)
                           * destFieldBox;
            return intersectionBox;
        }();

        // Per-direction centering is fixed over the box; parity varies per fine index. The full
        // stencil for any (centering, parity) combination is a distinct compile-time array; the
        // 2-bits-per-direction index built below (bit0 = parity, bit1 = dual) selects the matching
        // gather at runtime via gatherers_.
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

    // 1-D weight row for direction d and the compile-time choice (dual<<1 | parity). Forwards the
    // consteval primitive directly; rows differ in length (copy=1, half-point=2, dual-σ±=3) and
    // tensorProduct sizes the multi-D stencil from them exactly — no padding.
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
            value += w.coef * sampleCoarse_(src, anchorLocal, w.indexes);
        return value;
    }

    // runtime centering+parity → compile-time stencil: one gather fn per combination, dispatched
    // by a plain index. Each keeps its natural length; nothing is padded to a common size.
    using Gatherer_t                 = double (*)(FieldT const&, Point_t const&);
    static constexpr auto gatherers_ = core::for_N_make_array<nCombined>(
        [](auto ic) { return Gatherer_t{&gatherStencil_<ic()>}; });

    static double sampleCoarse_(FieldT const& src, Point_t const& anchorLocal,
                                Point_t const& offset)
    {
        if constexpr (dimension == 1)
            return src(anchorLocal[0] + offset[0]);
        else if constexpr (dimension == 2)
            return src(anchorLocal[0] + offset[0], anchorLocal[1] + offset[1]);
        else
            return src(anchorLocal[0] + offset[0], anchorLocal[1] + offset[1],
                       anchorLocal[2] + offset[2]);
    }

    // preserve the legacy NaN-guard: only fill fine indices not already set
    static void assignFine_(FieldT& dst, Point_t const& fineLocal, double const value)
    {
        if (auto& dst_val = dst(fineLocal.toArray()); std::isnan(dst_val))
            dst_val = value;
    }
};


// ---- factory (declared in field_refiner_kernel.hpp) ---------------------------------------------

template<typename GridLayoutT, typename FieldT>
std::shared_ptr<IFieldRefineKernel<GridLayoutT, FieldT>> makeRefineKernel(int const order)
{
    if (order != 2)
        throw std::runtime_error("makeRefineKernel: order must be 2 (Linear)");
    return std::make_shared<CompositeFieldRefiner<GridLayoutT, FieldT, 2>>();
}


} // namespace PHARE::amr


#endif
