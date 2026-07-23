#ifndef PHARE_COMPOSITE_FIELD_REFINER_HPP
#define PHARE_COMPOSITE_FIELD_REFINER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "field_refiner_kernel.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/IntVector.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>


namespace PHARE::amr
{

/** @brief Tag for the unlimited (smooth) prolongation path. */
struct NoLimiter
{
};


/**
 * @brief Runtime field-refinement kernel built from the two 1-D primitives tensor-producted.
 *
 * For ratio 2, each fine index maps to a coarse anchor I = floor(f/2) (toCoarseIndex) and a parity
 * p = f − 2I ∈ {0,1} per direction. Per direction the 1-D weight row is chosen by centering+parity:
 *   - primal, p=0 : exact copy (coincident node)            → {offset 0, 1.0}
 *   - primal, p=1 : half-point midpoint (Primitive A.2)     → directionalInterp<dir,PrimalToDual,order>
 *   - dual,   p   : ±¼ child ladder (Primitive B)           → directionalProlongation<dir, σ=2p−1, order>
 * The multi-D stencil is the outer product of the rows (tensorProductRuntime). Centering is constant
 * over a refineBox call, so the 2^dim parity-combo stencils are built once per box, then looked up
 * per fine index. Offsets are relative to the anchor I; values are gathered from the coarse field.
 */
template<typename GridLayoutT, typename FieldT, std::size_t order, bool sharedFacesOnly = false>
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

        // For a face field (B): the single primal direction is the face normal. Shared faces are
        // the even-normal-parity ones (coincident with a coarse face); the odd-normal (interior)
        // faces are owned by the Tóth-Roe postprocess, which overwrites them and never reads them,
        // so we skip them here. Tangential directions are dual ⇒ their children are always filled.
        // A B component is primal in at most one direction (its face normal). When it has one
        // (in-plane Bx/By/Bz), the odd-normal interior faces are owned by the Tóth-Roe postprocess
        // and skipped below. When it has none (out-of-plane / reduced-dimension component: By,Bz in
        // 1D, Bz in 2D), there is no normal face and no ∇·B interior to protect, so it prolongs like
        // any same-centered (all-dual) quantity — every fine face is filled.
        std::size_t normalDir = 0;
        bool hasNormal        = false;
        if constexpr (sharedFacesOnly)
        {
            int primalCount = 0;
            for (std::size_t d = 0; d < dimension; ++d)
                if (centering[d] == core::QtyCentering::primal)
                {
                    normalDir = d;
                    ++primalCount;
                }
            if (primalCount > 1)
                throw std::runtime_error(
                    "magnetic shared-face refiner expects at most one primal (normal) direction");
            hasNormal = (primalCount == 1);
        }

        // Centering is fixed over the box, so the rows are data-independent.
        // Build the 2^dim parity-combo stencils once, then look up per fine index.
        std::array<std::vector<WeightPoint_t>, (1u << dimension)> combos;
        for (std::size_t combo = 0; combo < combos.size(); ++combo)
            combos[combo] = makeStencil_(centering, combo);

        for (auto const fineIndex : phare_box_from<dimension>(intersectionBox))
        {
            auto const anchor = toCoarseIndex<dimension>(fineIndex);

            std::size_t combo = 0;
            for (std::size_t d = 0; d < dimension; ++d)
                combo |= static_cast<std::size_t>(fineIndex[d] - 2 * anchor[d]) << d;

            if constexpr (sharedFacesOnly)
                if (hasNormal && ((combo >> normalDir) & 1u) != 0u)
                    continue; // interior-normal face: left to Tóth-Roe postprocess

            auto const anchorLocal = AMRToLocal(anchor, sourceFieldBox);
            auto const fineLocal   = AMRToLocal(fineIndex, destFieldBox);

            double value = 0.;
            for (auto const& w : combos[combo])
                value += w.coef * sampleCoarse_(sourceField, anchorLocal, w.indexes);

            assignFine_(destinationField, fineLocal, value);
        }
    }

    // order 2 reads ±1 coarse cell (max |offset| over both 1-D primitives).
    int coarseStencilWidth() const override { return order / 2; }

private:
    // 1-D weight row for direction d, given its centering and parity (0/1), as a runtime vector.
    template<std::size_t d>
    static std::vector<WeightPoint_t> oneDRow_(core::QtyCentering c, std::size_t parity)
    {
        std::vector<WeightPoint_t> row;
        if (c == core::QtyCentering::primal)
        {
            if (parity == 0)
                row.push_back({Point_t{}, 1.0}); // coincident node: exact copy
            else
                for (auto const& w :
                     GridLayoutImpl::template directionalInterp<d, GridLayoutImpl::InterpDir::PrimalToDual, order>())
                    row.push_back(w);
        }
        else // dual: σ = 2·parity − 1
        {
            if (parity == 0)
                for (auto const& w : GridLayoutImpl::template directionalProlongation<d, -1, order>())
                    row.push_back(w);
            else
                for (auto const& w : GridLayoutImpl::template directionalProlongation<d, +1, order>())
                    row.push_back(w);
        }
        return row;
    }

    static std::vector<WeightPoint_t>
    makeStencil_(std::array<core::QtyCentering, dimension> const& centering, std::size_t combo)
    {
        auto const parity = [combo](std::size_t d) { return (combo >> d) & 1u; };

        if constexpr (dimension == 1)
            return oneDRow_<0>(centering[0], parity(0));
        else if constexpr (dimension == 2)
            return GridLayoutImpl::template tensorProductRuntime<0, 1>(
                oneDRow_<0>(centering[0], parity(0)), oneDRow_<1>(centering[1], parity(1)));
        else
            return GridLayoutImpl::template tensorProductRuntime<0, 1, 2>(
                oneDRow_<0>(centering[0], parity(0)), oneDRow_<1>(centering[1], parity(1)),
                oneDRow_<2>(centering[2], parity(2)));
    }

    static double sampleCoarse_(FieldT const& src, Point_t const& anchorLocal, Point_t const& offset)
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
    static void assignFine_(FieldT& dst, Point_t const& fineLocal, double value)
    {
        if constexpr (dimension == 1)
        {
            if (std::isnan(dst(fineLocal[0])))
                dst(fineLocal[0]) = value;
        }
        else if constexpr (dimension == 2)
        {
            if (std::isnan(dst(fineLocal[0], fineLocal[1])))
                dst(fineLocal[0], fineLocal[1]) = value;
        }
        else
        {
            if (std::isnan(dst(fineLocal[0], fineLocal[1], fineLocal[2])))
                dst(fineLocal[0], fineLocal[1], fineLocal[2]) = value;
        }
    }
};


// ---- factory (declared in field_refiner_kernel.hpp) ---------------------------------------------

template<typename GridLayoutT, typename FieldT>
std::shared_ptr<IFieldRefineKernel<GridLayoutT, FieldT>>
makeRefineKernel(int order)
{
    switch (order)
    {
        case 2: return std::make_shared<CompositeFieldRefiner<GridLayoutT, FieldT, 2>>();
        default: throw std::runtime_error("makeRefineKernel: order must be 2 (Linear)");
    }
}


} // namespace PHARE::amr


#endif
