#ifndef PHARE_COMPOSITE_FIELD_REFINER_HPP
#define PHARE_COMPOSITE_FIELD_REFINER_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"

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
namespace detail
{
    /**
     * @brief Where a fine index sits along one direction relative to its coarse anchor
     * I = floor(f/2), for refinement ratio 2.
     *
     * Along a primal direction the even fine node coincides with coarse node I and the odd one
     * lies halfway between nodes I and I+1. Along a dual direction coarse cell I splits into a
     * left (even) and a right (odd) fine cell.
     */
    enum class FineSite { OnCoarseNode, BetweenCoarseNodes, LeftChild, RightChild };

    NO_DISCARD constexpr FineSite fineSite(core::QtyCentering const centering, int const fineIndex)
    {
        bool const odd = isOddIndex(fineIndex);

        if (centering == core::QtyCentering::primal)
            return odd ? FineSite::BetweenCoarseNodes : FineSite::OnCoarseNode;

        return odd ? FineSite::RightChild : FineSite::LeftChild;
    }


    /**
     * @brief The 1-D prolongation weights of one FineSite: along that direction the fine value is
     * Σ_i weight[i] · coarse(I + offset[i]).
     */
    template<std::size_t order>
    struct ProlongationRow
    {
        // the widest row, a dual child, reads offsets -order/2 .. +order/2
        static constexpr std::size_t capacity = order + 1;

        std::array<int, capacity> offset{};
        std::array<double, capacity> weight{};
        std::size_t size = 0;
    };

    template<std::size_t order, std::size_t n>
    NO_DISCARD consteval ProlongationRow<order>
    toProlongationRow(std::array<core::WeightPoint<1>, n> const& weightPoints)
    {
        static_assert(n <= ProlongationRow<order>::capacity);

        ProlongationRow<order> row;
        for (auto const& wp : weightPoints)
        {
            row.offset[row.size] = wp.indexes[0];
            row.weight[row.size] = wp.coef;
            ++row.size;
        }
        return row;
    }


    /**
     * @brief The composite refinement fill: every fine index of a region from the coarse field.
     *
     * The stencil is separable. Along each direction the fine index's FineSite selects a 1-D row,
     *   - OnCoarseNode         : copy of node I
     *   - BetweenCoarseNodes   : midpoint of nodes I and I+1 (directionalInterp, PrimalToDual)
     *   - LeftChild, RightChild: the ±¼ child of cell I, which means back to it
     *                            (directionalProlongation)
     * and the fine value is the tensor product of those rows applied to the coarse field.
     *
     * Stateless: the kernels below are thin IFieldRefineKernel adapters over it that differ only
     * in the region they hand it, so a kernel costs one virtual call per overlap box.
     */
    template<typename GridLayoutT, typename FieldT, std::size_t order>
    struct CompositeFill
    {
        static constexpr std::size_t dimension = GridLayoutT::dimension;

        static_assert(order == 2, "composite refinement is order 2 (Linear)");

        // order 2 reads ±1 coarse cell, the midpoint and the child rows alike
        static constexpr int coarseStencilWidth() { return order / 2; }

        /**
         * @brief Fill every fine index of `region` from the coarse field.
         *
         * `region` is the caller's to choose: it must lie inside the destination allocation, and
         * only its fine indices are written.
         */
        static void fill(FieldT const& sourceField, FieldT& destinationField,
                         SAMRAI::hier::Box const& region,
                         std::array<core::QtyCentering, dimension> const& centering,
                         SAMRAI::hier::Box const& destFieldBox,
                         SAMRAI::hier::Box const& sourceFieldBox,
                         SAMRAI::hier::IntVector const& ratio)
        {
            for (std::size_t d = 0; d < dimension; ++d)
                if (ratio(d) != 2)
                    throw std::runtime_error(
                        "composite field refinement supports refinement ratio 2 only");

            for (auto const fineIndex : phare_box_from<dimension>(region))
            {
                std::array<Row const*, dimension> rows;
                for (std::size_t d = 0; d < dimension; ++d)
                    rows[d] = &row_(fineSite(centering[d], fineIndex[d]));

                auto const anchor = AMRToLocal(toCoarseIndex<dimension>(fineIndex), sourceFieldBox);
                auto const value  = gather_(sourceField, anchor, rows);

                assignIfUnset_(destinationField, AMRToLocal(fineIndex, destFieldBox), value);
            }
        }

    private:
        using Point_t = core::Point<int, dimension>;
        using Row     = ProlongationRow<order>;

        static constexpr Row onCoarseNode_
            = toProlongationRow<order>(std::array{core::WeightPoint<1>{core::Point<int, 1>{}, 1.}});
        static constexpr Row betweenCoarseNodes_ = toProlongationRow<order>(
            core::directionalInterp<1, 0, core::InterpDir::PrimalToDual, order>());
        static constexpr Row leftChild_
            = toProlongationRow<order>(core::directionalProlongation<1, 0, -1, order>());
        static constexpr Row rightChild_
            = toProlongationRow<order>(core::directionalProlongation<1, 0, +1, order>());

        static Row const& row_(FineSite const site)
        {
            switch (site)
            {
                case FineSite::OnCoarseNode: return onCoarseNode_;
                case FineSite::BetweenCoarseNodes: return betweenCoarseNodes_;
                case FineSite::LeftChild: return leftChild_;
                case FineSite::RightChild: return rightChild_;
            }
            throw std::logic_error("unknown FineSite");
        }

        // the tensor product of the rows, applied around the coarse anchor
        static double gather_(FieldT const& src, Point_t const& anchor,
                              std::array<Row const*, dimension> const& rows)
        {
            auto const& x = *rows[0];
            double value  = 0.;

            if constexpr (dimension == 1)
            {
                for (std::size_t i = 0; i < x.size; ++i)
                    value += x.weight[i] * src(anchor[0] + x.offset[i]);
            }
            else if constexpr (dimension == 2)
            {
                auto const& y = *rows[1];
                for (std::size_t i = 0; i < x.size; ++i)
                    for (std::size_t j = 0; j < y.size; ++j)
                        value += (x.weight[i] * y.weight[j])
                                 * src(anchor[0] + x.offset[i], anchor[1] + y.offset[j]);
            }
            else
            {
                auto const& y = *rows[1];
                auto const& z = *rows[2];
                for (std::size_t i = 0; i < x.size; ++i)
                    for (std::size_t j = 0; j < y.size; ++j)
                        for (std::size_t k = 0; k < z.size; ++k)
                            value += (x.weight[i] * y.weight[j] * z.weight[k])
                                     * src(anchor[0] + x.offset[i], anchor[1] + y.offset[j],
                                           anchor[2] + z.offset[k]);
            }

            return value;
        }

        // Writes only fine indices still NaN, so a value already in place (a same-level copy, or
        // an earlier overlapping fill) is never overwritten.
        static void assignIfUnset_(FieldT& dst, Point_t const& fineLocal, double const value)
        {
            if (auto& dst_val = dst(fineLocal); std::isnan(dst_val))
                dst_val = value;
        }
    };

} // namespace detail


/**
 * @brief Runtime field-refinement kernel: the composite fill over exactly the overlap SAMRAI hands
 * it.
 *
 * Every quantity but B refines this way. B needs its gather widened to whole coarse cells, which
 * is its own kernel (magnetic_composite_refiner.hpp).
 */
template<typename GridLayoutT, typename FieldT, std::size_t order>
class CompositeFieldRefiner final : public IFieldRefineKernel<GridLayoutT, FieldT>
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
        Fill::fill(sourceField, destinationField, intersectionBox, centering, destFieldBox,
                   sourceFieldBox, ratio);
    }

    int coarseStencilWidth() const override { return Fill::coarseStencilWidth(); }
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
