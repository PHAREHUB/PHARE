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

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>


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
     * @brief Call `fn` with `site` as a compile-time constant, so a stencil built from it stays a
     * constexpr table.
     */
    template<typename Fn>
    decltype(auto) withSite(FineSite const site, Fn&& fn)
    {
        using enum FineSite;

        switch (site)
        {
            case OnCoarseNode: return fn(std::integral_constant<FineSite, OnCoarseNode>{});
            case BetweenCoarseNodes:
                return fn(std::integral_constant<FineSite, BetweenCoarseNodes>{});
            case LeftChild: return fn(std::integral_constant<FineSite, LeftChild>{});
            case RightChild: return fn(std::integral_constant<FineSite, RightChild>{});
        }
        throw std::logic_error("unknown FineSite");
    }


    /**
     * @brief The composite refinement fill: every fine index of a region from the coarse field.
     *
     * The stencil is separable. Along each direction the fine index's FineSite selects a 1-D row,
     *   - OnCoarseNode         : copy of node I
     *   - BetweenCoarseNodes   : midpoint of nodes I and I+1 (directionalInterp, PrimalToDual)
     *   - LeftChild, RightChild: the ±¼ child of cell I, which means back to it
     *                            (directionalProlongation)
     * and the fine value is the tensor product of those rows applied to the coarse field. Every
     * combination of sites is a compile-time stencil, selected per fine index.
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
                std::array<FineSite, dimension> fineSites;
                for (std::size_t d = 0; d < dimension; ++d)
                    fineSites[d] = fineSite(centering[d], fineIndex[d]);

                auto const anchor = AMRToLocal(toCoarseIndex<dimension>(fineIndex), sourceFieldBox);
                auto const value  = gather_<>(sourceField, anchor, fineSites);

                assignIfUnset_(destinationField, AMRToLocal(fineIndex, destFieldBox), value);
            }
        }

    private:
        using Point_t = core::Point<int, dimension>;

        // the 1-D row of direction `dir` for a fine index at `site`
        template<std::size_t dir, FineSite site>
        static consteval auto row_()
        {
            if constexpr (site == FineSite::OnCoarseNode)
                return std::array{core::WeightPoint<dimension>{Point_t{}, 1.}};
            else if constexpr (site == FineSite::BetweenCoarseNodes)
                return core::directionalInterp<dimension, dir, core::InterpDir::PrimalToDual,
                                               order>();
            else if constexpr (site == FineSite::LeftChild)
                return core::directionalProlongation<dimension, dir, -1, order>();
            else
                return core::directionalProlongation<dimension, dir, +1, order>();
        }

        // Recursive: each call moves the next direction's runtime fine site into the compile-time
        // `site` pack. Once every direction is in it, the stencil is the constexpr tensor product
        // of their rows, summed around the coarse anchor.
        template<FineSite... site>
        static double gather_(FieldT const& src, Point_t const& anchor,
                              std::array<FineSite, dimension> const& fineSites)
        {
            if constexpr (sizeof...(site) == dimension)
            {
                return core::for_N_fold<dimension>([&](auto... dir) {
                    constexpr auto stencil = core::tensorProduct<decltype(dir)::value...>(
                        row_<decltype(dir)::value, site>()...);

                    return core::sum_from(
                        stencil, [&](auto const& w) { return w.coef * src(anchor + w.indexes); });
                });
            }
            else
            {
                return withSite(fineSites[sizeof...(site)], [&](auto next) {
                    return gather_<site..., decltype(next)::value>(src, anchor, fineSites);
                });
            }
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
