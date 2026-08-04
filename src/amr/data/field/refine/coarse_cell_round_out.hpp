#ifndef PHARE_AMR_COARSE_CELL_ROUND_OUT_HPP
#define PHARE_AMR_COARSE_CELL_ROUND_OUT_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"

#include <SAMRAI/hier/Box.h>

#include <array>
#include <cstddef>
#include <string>

namespace PHARE::amr
{
/**
 * @brief The whole-coarse-cell invariant of magnetic prolongation, and the helpers enforcing it.
 *
 * Magnetic prolongation (ratio 2) is split in two. Faces with an even index in their component's
 * normal direction coincide with a coarse face and are gathered from the coarse level; the
 * odd-normal (interior) faces are then reconstructed div-free by a magnetic patch strategy's
 * postprocessRefine (see MagneticPatchStrategyBase::reconstructionRegion). That reconstruction
 * reaches exactly one coarse cell: for an interior face of coarse cell C it reads only faces
 * bounding C, and every one of them is a shared face the gather just filled.
 *
 * Hence the invariant: **the region a magnetic prolongation runs over must be a union of whole
 * coarse cells**. Then every postprocess input is in the region and was written. SAMRAI's fill
 * boxes do not satisfy it on their own: a recursive schedule fills a coarse-interpolation
 * temporary plus a ring of d_max_stencil_width = 1 cell, and one cell is always *half* a coarse
 * cell.
 *
 * A cell box is a union of whole coarse cells in a direction iff its lower index is even and its
 * upper index is odd. The corresponding field box, per direction, is
 *   - dual  : rows  2I .. 2J+1  (lower even, upper odd)
 *   - primal: nodes 2I .. 2J+2  (lower even, upper even)
 *
 * Parity is tested with `% 2 != 0`, not `== 1`: AMR indices are signed and routinely negative, and
 * `i % 2` is -1 there.
 */

NO_DISCARD constexpr bool isOddIndex(int const i)
{
    return i % 2 != 0;
}

NO_DISCARD constexpr int roundDownToEven(int const i)
{
    return isOddIndex(i) ? i - 1 : i;
}

NO_DISCARD constexpr int roundUpToOddIndex(int const i)
{
    return isOddIndex(i) ? i : i + 1;
}

NO_DISCARD constexpr int roundUpToEvenIndex(int const i)
{
    return isOddIndex(i) ? i + 1 : i;
}


//! smallest union of whole coarse cells containing the given cell box
template<std::size_t dim>
NO_DISCARD SAMRAI::hier::Box roundCellBoxOutToCoarseCells(SAMRAI::hier::Box box)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        auto const dir = static_cast<SAMRAI::hier::Box::dir_t>(d);
        box.setLower(dir, roundDownToEven(box.lower(dir)));
        box.setUpper(dir, roundUpToOddIndex(box.upper(dir)));
    }
    return box;
}


//! the same round-out expressed on a field box, i.e. the field box of the rounded-out cell box
template<std::size_t dim>
NO_DISCARD SAMRAI::hier::Box
roundFieldBoxOutToCoarseCells(SAMRAI::hier::Box box,
                              std::array<core::QtyCentering, dim> const& centering)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        auto const dir = static_cast<SAMRAI::hier::Box::dir_t>(d);
        box.setLower(dir, roundDownToEven(box.lower(dir)));
        box.setUpper(dir, centering[d] == core::QtyCentering::primal
                              ? roundUpToEvenIndex(box.upper(dir))
                              : roundUpToOddIndex(box.upper(dir)));
    }
    return box;
}


/**
 * @brief Does this cell box consist of whole coarse cells, i.e. does it satisfy the invariant?
 *
 * Checking a region once is equivalent to checking every index in it: if the region is a union of
 * whole coarse cells then the coarse cell of every fine face it contains is inside it too, which is
 * exactly what each Tóth-Roe reconstruction needs.
 */
template<std::size_t dim>
NO_DISCARD bool isWholeCoarseCells(SAMRAI::hier::Box const& box)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        auto const dir = static_cast<SAMRAI::hier::Box::dir_t>(d);
        if (isOddIndex(box.lower(dir)) or not isOddIndex(box.upper(dir)))
            return false;
    }
    return true;
}


//! for diagnostics: "[(lo0,lo1),(up0,up1)]"
NO_DISCARD inline std::string to_str(SAMRAI::hier::Box const& box)
{
    auto const dim = box.getDim().getValue();
    std::string s{"[("};
    for (SAMRAI::hier::Box::dir_t d = 0; d < dim; ++d)
        s += std::to_string(box.lower(d)) + (d + 1 < dim ? "," : "),(");
    for (SAMRAI::hier::Box::dir_t d = 0; d < dim; ++d)
        s += std::to_string(box.upper(d)) + (d + 1 < dim ? "," : ")]");
    return s;
}

} // namespace PHARE::amr

#endif // PHARE_AMR_COARSE_CELL_ROUND_OUT_HPP
