#ifndef PHARE_AMR_ADPT_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_ADPT_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/constants.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/refine/coarse_cell_round_out.hpp"

#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

/**
 * @brief Stage 2 of the Balsara ADPT divergence-free magnetic prolongation: the cross-component
 *        divB touch-up. Order-independent.
 *
 * Stage 1 (the fill-all composite kernel, magnetic_composite_refiner.hpp) fills every fine face of
 * each B component from that component's own coarse faces, and makes no divB claim. This stage
 * adds a closed-form min-norm correction to the 2d interior faces of each coarse zone, equalizing
 * the divergences of the 2^d fine cells the zone splits into -- its subzones. They all become the
 * zone's transported divergence q0, so a discretely div-free coarse field (q0 = 0) prolongs
 * div-free to roundoff whatever the stage-1 order.
 *
 * Balsara, D. S., Samantaray, S. & Subramanian, S. (2023), "Efficient WENO-Based Prolongation
 * Strategies for Divergence-Preserving Vector Fields", Commun. Appl. Math. Comput.,
 * doi:10.1007/s42967-021-00182-x. The correction of an interior face is a weighted sum of the
 * zone's subzone divergence deficits — own subzone weighted highest, each farther
 * (Hamming-distance) neighbour less. The per-dimension weights are the flux-variable ones,
 * 1/2 (1D), [3,1]/8 (2D) and [7,2,2,1]/24 (3D): the density weights ([·]/16, [·]/48) rescaled by
 * the refinement ratio 2.
 *
 * The correction is added to the stage-1 face value rather than replacing it, so a higher-order
 * stage-1 interior fill survives the touch-up.
 */
template<typename TensorFieldDataT>
class ADPTMagneticRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    using Geometry        = typename TensorFieldDataT::Geometry;
    using gridlayout_type = typename TensorFieldDataT::gridlayout_type;

    static constexpr std::size_t N         = TensorFieldDataT::N;
    static constexpr std::size_t dimension = TensorFieldDataT::dimension;

    using CellKey = std::array<int, dimension>;

    /**
     * @brief The stage-1 divergence snapshot of one touch-up pass: one slot per fine cell of the
     * reconstruction region, NaN meaning "not computed yet".
     *
     * Every correction has to read stage-1 divergences -- recomputing from the live field once a
     * sibling interior face has been written couples the two faces and loses divB exactness -- so
     * each subzone divergence is computed once per pass and reread afterwards. The region is a
     * dense box of fine cells, so the slots are a flat array indexed by the cell's offset in it,
     * and the array is kept between passes: only a region bigger than every earlier one allocates.
     */
    class DivScratch
    {
    public:
        //! sizes the snapshot to `region`'s fine cells, in the local indices the corrections
        //! index by, and marks every slot not-computed
        void reset(SAMRAI::hier::Box const& region, gridlayout_type const& layout)
        {
            auto const local = layout.AMRToLocal(phare_box_from<dimension>(region));

            lower_ = local.lower;
            shape_ = local.shape();
            values_.assign(local.size(), unset_);
        }

        NO_DISCARD double& operator()(auto const... cells)
        {
            static_assert(sizeof...(cells) == dimension);

            return values_[offset_(CellKey{cells...})];
        }

    private:
        static constexpr double unset_ = std::numeric_limits<double>::quiet_NaN();

        // row-major offset of a fine cell in the region. The corrections only ever reach cells of
        // the coarse cell they correct an interior face of, so a region of whole coarse cells
        // (reconstructionRegion) contains every cell they index -- the assert catches a caller
        // that reset the snapshot to something else.
        NO_DISCARD std::size_t offset_(CellKey const& cell) const
        {
            std::size_t offset = 0;

            for (std::size_t d = 0; d < dimension; ++d)
            {
                assert(cell[d] >= lower_[d] and cell[d] - lower_[d] < shape_[d]);

                offset = offset * static_cast<std::size_t>(shape_[d])
                         + static_cast<std::size_t>(cell[d] - lower_[d]);
            }

            return offset;
        }

        core::Point<int, dimension> lower_;
        core::Point<int, dimension> shape_;
        std::vector<double> values_;
    };

    ADPTMagneticRefinePatchStrategy()
        : b_id_{-1}
    {
    }

    void assertIDsSet() const
    {
        if (b_id_ < 0)
            throw std::runtime_error(
                "ADPTMagneticRefinePatchStrategy: registerIDs was not called before use");
    }

    void registerIDs(int const b_id) { b_id_ = b_id; }

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch&, double const,
                                       SAMRAI::hier::IntVector const&) override
    {
    }

    void preprocessRefine(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                          SAMRAI::hier::Box const&, SAMRAI::hier::IntVector const&) override
    {
    }

    // Always 0: the touch-up reads only fine faces the stage-1 gather already wrote, no coarse
    // data. SAMRAI provisions max(this, every registered refine operator), so 0 narrows nothing.
    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getZero(dim);
    }

    /**
     * @brief The region the touch-up over fill box `fine_box` must run over: `fine_box` rounded
     * out to whole coarse cells, clipped to the allocation.
     *
     * The touch-up reaches exactly one coarse cell, so it is self-consistent only over a union of
     * whole coarse cells (coarse_cell_round_out.hpp). SAMRAI fill boxes are not: a recursive
     * schedule's coarse-interpolation temporary plus its one-cell ring always cuts a coarse cell
     * in half, leaving the far shared faces unwritten.
     *
     * Round-out grows each edge by one index, so the region stays inside the allocation iff
     * field_ghost_width >= fill_ring + 1. That holds everywhere we support, but with zero slack:
     * the ring is 1 and field_ghost_width (ghost_width_calculator.hpp) bottoms out at exactly 2.
     * If the clip ever bit it would halve a coarse cell, and a half-covered cell has no well-posed
     * reconstruction — some of its inputs are faces nothing ever wrote — hence the throw.
     */
    static SAMRAI::hier::Box reconstructionRegion(SAMRAI::hier::Box const& fine_box,
                                                  SAMRAI::hier::Box const& ghost_box)
    {
        auto const region = roundCellBoxOutToCoarseCells<dimension>(fine_box) * ghost_box;

        if (!isWholeCoarseCells<dimension>(region))
            throw std::runtime_error(
                "magnetic prolongation region is not a union of whole coarse cells: fill box "
                + to_str(fine_box) + ", rounded out and clipped to the allocated "
                + to_str(ghost_box) + ", gives " + to_str(region)
                + ", which half-covers a coarse cell. The field ghost width must be at least the "
                  "fill ring plus one.");

        return region;
    }


    // Run over the whole-coarse-cell-rounded region so every subzoneDiv_ read lands on a fine
    // face the stage-1 refinement actually wrote.
    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
        assertIDsSet();

        auto& fields = TensorFieldDataT::getFields(fine, b_id_);

        auto const layout = PHARE::amr::layoutFromPatch<gridlayout_type>(fine);
        auto const region = reconstructionRegion(fine_box, fine.getPatchData(b_id_)->getGhostBox());

        touchUpInteriorFaces(fields, layout, region, scratch_);
    }


    /**
     * @brief Add the divergence-equalizing correction over `region`, which must already be
     * whole-coarse-cell rounded (reconstructionRegion): every read here lands on a fine face of
     * the coarse cell being corrected, and those exist only if the cell is wholly inside.
     */
    static void touchUpInteriorFaces(auto& fields, gridlayout_type const& layout,
                                     SAMRAI::hier::Box const& region, DivScratch& scratch)
    {
        auto& [bx, by, bz] = fields;

        auto const regionLayout = Geometry::layoutFromBox(region, layout);

        auto const fine_field_box = core::for_N_make_array<N>([&](auto i) {
            using PhysicalQuantity = std::decay_t<decltype(fields[i].physicalQuantity())>;

            return FieldGeometry<gridlayout_type, PhysicalQuantity>::toFieldBox(
                region, fields[i].physicalQuantity(), regionLayout);
        });

        // One stage-1 snapshot per pass (DivScratch): every correction must read stage-1
        // divergences, and the corrections write the faces those divergences are computed from.
        scratch.reset(region, layout);

        if constexpr (dimension == 1)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
                correctBx1d(scratch, bx, layout, i);
        }
        else if constexpr (dimension == 2)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
                correctBx2d(scratch, bx, by, layout, i);

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirY]))
                correctBy2d(scratch, bx, by, layout, i);
        }
        else if constexpr (dimension == 3)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
                correctBx3d(scratch, bx, by, bz, layout, i);

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirY]))
                correctBy3d(scratch, bx, by, bz, layout, i);

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirZ]))
                correctBz3d(scratch, bx, by, bz, layout, i);
        }
    }


    static auto isNewFineFace(auto const& amrIdx, auto const dir)
    {
        // amr index can be negative so test != 0 (odd) rather than == 1
        return amrIdx[dir] % 2 != 0;
    }


    // ---- 1D ------------------------------------------------------------------------------------
    // Only Bx has an x-normal; the single interior face's min-norm correction is δ = pair/2 (flux).
    // On div-free (Bx const) stage-1 data pair = 0, so this is a no-op.
    static void correctBx1d(auto& scratch, auto& bx, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirX))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];

        double const pair = subzoneDiv1d_(scratch, bx, ix) - subzoneDiv1d_(scratch, bx, ix - 1);
        bx(ix) += 0.5 * pair;
    }


    // ---- 2D ------------------------------------------------------------------------------------
    // Interior Bx face bx(ix,iy) separates fine cells (ix-1,iy) [left] and (ix,iy) [right].
    // ξ = [ 3·pair(own row) + 1·pair(sibling row) ] / 8   ([3,1] weights, flux denominator 8)
    // where pair(cy) = d(right cell) − d(left cell) = deficit difference across the face.
    static void correctBx2d(auto& scratch, auto& bx, auto& by, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirX))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];
        int const iy   = loc[dirY];

        int const cxL = ix - 1;
        int const cxR = ix;
        int const cy0 = iy;                                   // own row
        int const cy1 = iy + ((idx[dirY] % 2 == 0) ? 1 : -1); // sibling row in the coarse cell

        auto const& D = layout.meshSize();
        auto pair     = [&](int cy) {
            return subzoneDiv2d_(scratch, bx, by, D, cxR, cy)
                   - subzoneDiv2d_(scratch, bx, by, D, cxL, cy);
        };

        bx(ix, iy) += D[dirX] * (3.0 * pair(cy0) + pair(cy1)) / 8.0;
    }

    // Interior By face by(ix,iy) separates fine cells (ix,iy-1) [below] and (ix,iy) [above].
    // η = [ 3·pair(own column) + 1·pair(sibling column) ] / 8
    // where pair(cx) = d(above cell) − d(below cell).
    static void correctBy2d(auto& scratch, auto& bx, auto& by, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirY))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];
        int const iy   = loc[dirY];

        int const cyB = iy - 1;
        int const cyA = iy;
        int const cx0 = ix;                                   // own column
        int const cx1 = ix + ((idx[dirX] % 2 == 0) ? 1 : -1); // sibling column

        auto const& D = layout.meshSize();
        auto pair     = [&](int cx) {
            return subzoneDiv2d_(scratch, bx, by, D, cx, cyA)
                   - subzoneDiv2d_(scratch, bx, by, D, cx, cyB);
        };

        by(ix, iy) += D[dirY] * (3.0 * pair(cx0) + pair(cx1)) / 8.0;
    }


    // ---- 3D ------------------------------------------------------------------------------------
    // Interior face of component c on its midplane, transverse quadrant t: weights [7,2,2,1]/24
    // (flux). pair(t) = d(high cell along c) − d(low cell along c), evaluated at transverse
    // quadrant t; own quadrant weight 7, single-flip (edge-adjacent) 2 each, double-flip
    // (diagonal) 1 — all positive.
    static void correctBx3d(auto& scratch, auto& bx, auto& by, auto& bz, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirX))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];
        int const iy   = loc[dirY];
        int const iz   = loc[dirZ];

        int const cxL = ix - 1;
        int const cxR = ix;
        int const sy  = iy + ((idx[dirY] % 2 == 0) ? 1 : -1);
        int const sz  = iz + ((idx[dirZ] % 2 == 0) ? 1 : -1);

        auto const& D = layout.meshSize();
        auto pair     = [&](int cy, int cz) {
            return subzoneDiv3d_(scratch, bx, by, bz, D, cxR, cy, cz)
                   - subzoneDiv3d_(scratch, bx, by, bz, D, cxL, cy, cz);
        };

        bx(ix, iy, iz)
            += D[dirX]
               * (7.0 * pair(iy, iz) + 2.0 * pair(sy, iz) + 2.0 * pair(iy, sz) + pair(sy, sz))
               / 24.0;
    }

    static void correctBy3d(auto& scratch, auto& bx, auto& by, auto& bz, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirY))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];
        int const iy   = loc[dirY];
        int const iz   = loc[dirZ];

        int const cyB = iy - 1;
        int const cyA = iy;
        int const sx  = ix + ((idx[dirX] % 2 == 0) ? 1 : -1);
        int const sz  = iz + ((idx[dirZ] % 2 == 0) ? 1 : -1);

        auto const& D = layout.meshSize();
        auto pair     = [&](int cx, int cz) {
            return subzoneDiv3d_(scratch, bx, by, bz, D, cx, cyA, cz)
                   - subzoneDiv3d_(scratch, bx, by, bz, D, cx, cyB, cz);
        };

        by(ix, iy, iz)
            += D[dirY]
               * (7.0 * pair(ix, iz) + 2.0 * pair(sx, iz) + 2.0 * pair(ix, sz) + pair(sx, sz))
               / 24.0;
    }

    static void correctBz3d(auto& scratch, auto& bx, auto& by, auto& bz, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirZ))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];
        int const iy   = loc[dirY];
        int const iz   = loc[dirZ];

        int const czB = iz - 1;
        int const czA = iz;
        int const sx  = ix + ((idx[dirX] % 2 == 0) ? 1 : -1);
        int const sy  = iy + ((idx[dirY] % 2 == 0) ? 1 : -1);

        auto const& D = layout.meshSize();
        auto pair     = [&](int cx, int cy) {
            return subzoneDiv3d_(scratch, bx, by, bz, D, cx, cy, czA)
                   - subzoneDiv3d_(scratch, bx, by, bz, D, cx, cy, czB);
        };

        bz(ix, iy, iz)
            += D[dirZ]
               * (7.0 * pair(ix, iy) + 2.0 * pair(sx, iy) + 2.0 * pair(ix, sy) + pair(sx, sy))
               / 24.0;
    }


private:
    // 1/D_c-weighted divergence of the fine cell at local index (cx[,cy[,cz]]): the sum over
    // directions of (high face − low face)/D_c. 1D keeps the raw difference — its weight cancels
    // against the prefactor. Held in the pass snapshot so a correction still reads the stage-1
    // value after sibling faces are written.
    static double subzoneDiv1d_(auto& scratch, auto& bx, int cx)
    {
        auto& d = scratch(cx);

        if (std::isnan(d))
            d = bx(cx + 1) - bx(cx);

        return d;
    }

    static double subzoneDiv2d_(auto& scratch, auto& bx, auto& by, auto const& D, int cx, int cy)
    {
        auto& d = scratch(cx, cy);

        if (std::isnan(d))
            d = (bx(cx + 1, cy) - bx(cx, cy)) / D[dirX] + (by(cx, cy + 1) - by(cx, cy)) / D[dirY];

        return d;
    }

    static double subzoneDiv3d_(auto& scratch, auto& bx, auto& by, auto& bz, auto const& D, int cx,
                                int cy, int cz)
    {
        auto& d = scratch(cx, cy, cz);

        if (std::isnan(d))
            d = (bx(cx + 1, cy, cz) - bx(cx, cy, cz)) / D[dirX]
                + (by(cx, cy + 1, cz) - by(cx, cy, cz)) / D[dirY]
                + (bz(cx, cy, cz + 1) - bz(cx, cy, cz)) / D[dirZ];

        return d;
    }

    int b_id_;
    DivScratch scratch_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_ADPT_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
