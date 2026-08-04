#ifndef PHARE_AMR_ADPT_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_ADPT_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/constants.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/refine/magnetic_patch_strategy_base.hpp"

#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <array>
#include <cmath>
#include <cassert>
#include <map>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

/**
 * @brief Stage-2 (order-INdependent) cross-component divB touch-up of the Balsara ADPT
 *        div-free magnetic prolongation.
 *
 * Paired with the fill-all composite kernel (stage 1, per-component, order-dialed): stage 1
 * fills EVERY fine face of each B component from its own coarse faces; this stage equalizes,
 * per coarse zone, the 2^d fine-subzone divergences of the filled data by adding the closed-form
 * min-norm correction to the 2d interior faces. When the coarse field is discretely div-free the
 * transported zone divergence q0 = 0, so every corrected fine subzone divergence becomes 0 to
 * roundoff (div-free prolongation). At order 2 this touch-up is algebraically identical to the
 * legacy Tóth-Roe postprocess (derivation §S6b), so ADPT order 2 == the legacy operator; at order
 * 4 the same touch-up is ADDED to the 4th-order (Cubic4) interior faces, unlocking genuine
 * 4th-order div-free B refinement (out of scope here: this branch caps refinement order at 2).
 *
 * Derivation:
 * ~/.claude/plans/ucaromel/higher-order-refinement/2026-07-12-balsara-divfree-prolongation/
 *             derivation.md (§S4 2D, §S5 3D, §S10 implementation entry).
 *
 * IMPLEMENTATION NOTES
 * --------------------
 * - dx = dy (= dz) is ASSUMED (no runtime check, TR precedent — the legacy Laplacian correction
 *   silently shares this assumption). The mesh spacing cancels in the equal-mesh closed forms;
 *   see the §S6b DECISION. The anisotropic-mesh generalisation stays in the derivation reference.
 * - The correction is ADDED (+=) to the existing (stage-1) face value, whereas TR fully recomputed
 *   the interior face (baseline + correction). At order 2 the stage-1 interior fill equals TR's
 *   ½(neighbour) baseline, so += and TR's = coincide; at order 4 the += preserves the Cubic4 face.
 * - subzoneDiv here is the raw fine-face-difference sum (a "flux" divergence, no 1/h factor — the
 *   dimensionally consistent quantity to add to a face value). Consequently the closed-form
 *   denominators are the FLUX-variable ones: 1/2 (1D), [3,1]/8 (2D), [7,2,2,1]/24 (3D). These are
 *   the density weights of §S4/§S5 ([·]/16, [·]/48) rescaled by the fixed 1/h = 2 (the derivation's
 *   "flux variables ⇒ /8 and /24"); numerically reproduces TR exactly at equal mesh.
 * - Order-independence / correctness: the correction of every interior face must be computed from
 *   the STAGE-1 (uncorrected) subzone divergences. A face-centric loop that recomputed divergences
 *   from the live (partially corrected) field would couple sibling interior faces and fail the
 *   divB-exactness gate. A per-zone snapshot is therefore load-bearing (NOT an optional
 *   optimisation): subzoneDiv_ memoises each fine cell's divergence in a small per-postprocess
 *   cache, populated at the cell's first reference — which, since the touch-up only ever writes
 *   interior faces and a cell's divergence is only ever requested by its own zone's corrections,
 *   captures the stage-1 value.
 * - Whole-coarse-cell round-out: every subzoneDiv_ read of an interior face's touch-up reaches only
 *   the fine faces of its own coarse cell — the same one-coarse-cell reach the legacy Tóth-Roe
 *   postprocess has (see MagneticPatchStrategyBase::reconstructionRegion). So the touch-up must run
 *   over the SAME whole-coarse-cell region as the legacy strategy, computed the same way: round the
 *   SAMRAI fill box out to whole coarse cells and clip to the allocated ghost box. Running it over
 *   the raw fill_box (as an earlier version of this file did) reopens the 3-level perimeter-NaN bug
 *   the round-out fixed: on a misaligned fill box, the far shared face of a cut coarse cell was
 *   never gathered, and correctBxNd would read it as leftover NaN scratch.
 */
template<typename ResMan, typename TensorFieldDataT>
class ADPTMagneticRefinePatchStrategy : public MagneticPatchStrategyBase
{
public:
    using Geometry        = typename TensorFieldDataT::Geometry;
    using gridlayout_type = typename TensorFieldDataT::gridlayout_type;

    static constexpr std::size_t N         = TensorFieldDataT::N;
    static constexpr std::size_t dimension = TensorFieldDataT::dimension;

    using CellKey  = std::array<int, dimension>;
    using DivCache = std::map<CellKey, double>;

    ADPTMagneticRefinePatchStrategy(ResMan& resourcesManager)
        : rm_{resourcesManager}
        , b_id_{-1}
    {
    }

    void assertIDsSet() const
    {
        assert(b_id_ >= 0 && "ADPTMagneticRefinePatchStrategy: IDs must be registered before use");
    }

    void registerIDs(int const b_id) override { b_id_ = b_id; }

    //! forwards to the shared MagneticPatchStrategyBase implementation with this strategy's own
    //! (compile-time) dimension; see the base for the invariant and the reasoning.
    static SAMRAI::hier::Box reconstructionRegion(SAMRAI::hier::Box const& fine_box,
                                                  SAMRAI::hier::Box const& ghost_box)
    {
        return MagneticPatchStrategyBase::reconstructionRegion<dimension>(fine_box, ghost_box);
    }


    // Stage-2 touch-up: same postprocessRefine slot as the legacy TR strategy, run over the same
    // whole-coarse-cell-rounded region (reconstructionRegion) so every subzoneDiv_ read lands on a
    // fine face the stage-1 gather actually wrote (see class notes).
    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
        assertIDsSet();

        auto& fields = TensorFieldDataT::getFields(fine, b_id_);

        auto const layout = PHARE::amr::layoutFromPatch<gridlayout_type>(fine);
        auto const region = reconstructionRegion(fine_box, fine.getPatchData(b_id_)->getGhostBox());

        touchUpInteriorFaces(fields, layout, region);
    }


    /**
     * @brief Stage 2 proper: add the divergence-equalizing correction over `region`.
     *
     * Split out of postprocessRefine so it can be driven directly from a test with a hand-built
     * field, exactly as the legacy strategy's reconstructInteriorFaces was. `region` must already
     * be whole-coarse-cell rounded (reconstructionRegion) — every read here lands on a fine face
     * of the coarse cell being corrected, and those exist only if the cell is wholly inside.
     */
    static void touchUpInteriorFaces(auto& fields, gridlayout_type const& layout,
                                     SAMRAI::hier::Box const& region)
    {
        auto& [bx, by, bz] = fields;

        auto const regionLayout = Geometry::layoutFromBox(region, layout);

        auto const fine_field_box = core::for_N_make_array<N>([&](auto i) {
            using PhysicalQuantity = std::decay_t<decltype(fields[i].physicalQuantity())>;

            return FieldGeometry<gridlayout_type, PhysicalQuantity>::toFieldBox(
                region, fields[i].physicalQuantity(), regionLayout);
        });

        // one stage-1 divergence snapshot for the whole patch pair (see class notes)
        DivCache cache;

        if constexpr (dimension == 1)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
                correctBx1d(cache, bx, layout, i);
        }
        else if constexpr (dimension == 2)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
                correctBx2d(cache, bx, by, layout, i);

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirY]))
                correctBy2d(cache, bx, by, layout, i);
        }
        else if constexpr (dimension == 3)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
                correctBx3d(cache, bx, by, bz, layout, i);

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirY]))
                correctBy3d(cache, bx, by, bz, layout, i);

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirZ]))
                correctBz3d(cache, bx, by, bz, layout, i);
        }
    }


    static auto isNewFineFace(auto const& amrIdx, auto const dir)
    {
        // amr index can be negative so test != 0 (odd) rather than == 1
        return amrIdx[dir] % 2 != 0;
    }


    // ---- 1D ------------------------------------------------------------------------------------
    // Only Bx has an x-normal; the single interior face's min-norm correction is δ = pair/2 (flux).
    // On div-free (Bx const) stage-1 data pair = 0, so this is a no-op — matching TR's 1D baseline.
    static void correctBx1d(auto& cache, auto& bx, auto const& layout,
                            core::Point<int, dimension> idx)
    {
        if (!isNewFineFace(idx, dirX))
            return;

        auto const loc = layout.AMRToLocal(idx);
        int const ix   = loc[dirX];

        double const pair = subzoneDiv1d_(cache, bx, ix) - subzoneDiv1d_(cache, bx, ix - 1);
        bx(ix) += 0.5 * pair;
    }


    // ---- 2D ------------------------------------------------------------------------------------
    // Interior Bx face bx(ix,iy) separates fine cells (ix-1,iy) [left] and (ix,iy) [right].
    // ξ = [ 3·pair(own row) + 1·pair(sibling row) ] / 8   (§S4 [3,1], flux denominator 8)
    // where pair(cy) = d(right cell) − d(left cell) = deficit difference across the face.
    static void correctBx2d(auto& cache, auto& bx, auto& by, auto const& layout,
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

        auto pair = [&](int cy) {
            return subzoneDiv2d_(cache, bx, by, cxR, cy) - subzoneDiv2d_(cache, bx, by, cxL, cy);
        };

        bx(ix, iy) += (3.0 * pair(cy0) + pair(cy1)) / 8.0;
    }

    // Interior By face by(ix,iy) separates fine cells (ix,iy-1) [below] and (ix,iy) [above].
    // η = [ 3·pair(own column) + 1·pair(sibling column) ] / 8
    // where pair(cx) = d(above cell) − d(below cell).
    static void correctBy2d(auto& cache, auto& bx, auto& by, auto const& layout,
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

        auto pair = [&](int cx) {
            return subzoneDiv2d_(cache, bx, by, cx, cyA) - subzoneDiv2d_(cache, bx, by, cx, cyB);
        };

        by(ix, iy) += (3.0 * pair(cx0) + pair(cx1)) / 8.0;
    }


    // ---- 3D ------------------------------------------------------------------------------------
    // Interior face of component c on its midplane, transverse quadrant t: weights [7,2,2,1]/24
    // (§S5, flux). pair(t) = d(high cell along c) − d(low cell along c), evaluated at transverse
    // quadrant t; own quadrant weight 7, single-flip (edge-adjacent) 2 each, double-flip
    // (diagonal) 1 — all positive.
    static void correctBx3d(auto& cache, auto& bx, auto& by, auto& bz, auto const& layout,
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

        auto pair = [&](int cy, int cz) {
            return subzoneDiv3d_(cache, bx, by, bz, cxR, cy, cz)
                   - subzoneDiv3d_(cache, bx, by, bz, cxL, cy, cz);
        };

        bx(ix, iy, iz)
            += (7.0 * pair(iy, iz) + 2.0 * pair(sy, iz) + 2.0 * pair(iy, sz) + pair(sy, sz)) / 24.0;
    }

    static void correctBy3d(auto& cache, auto& bx, auto& by, auto& bz, auto const& layout,
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

        auto pair = [&](int cx, int cz) {
            return subzoneDiv3d_(cache, bx, by, bz, cx, cyA, cz)
                   - subzoneDiv3d_(cache, bx, by, bz, cx, cyB, cz);
        };

        by(ix, iy, iz)
            += (7.0 * pair(ix, iz) + 2.0 * pair(sx, iz) + 2.0 * pair(ix, sz) + pair(sx, sz)) / 24.0;
    }

    static void correctBz3d(auto& cache, auto& bx, auto& by, auto& bz, auto const& layout,
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

        auto pair = [&](int cx, int cy) {
            return subzoneDiv3d_(cache, bx, by, bz, cx, cy, czA)
                   - subzoneDiv3d_(cache, bx, by, bz, cx, cy, czB);
        };

        bz(ix, iy, iz)
            += (7.0 * pair(ix, iy) + 2.0 * pair(sx, iy) + 2.0 * pair(ix, sy) + pair(sx, sy)) / 24.0;
    }


private:
    // Raw ("flux") divergence of the fine cell at local index (cx[,cy[,cz]]): the sum over
    // directions of the (high face − low face) difference (dx = dy = dz assumed). Memoised so every
    // interior face's correction reads the STAGE-1 value even after sibling faces are written.
    static double subzoneDiv1d_(auto& cache, auto& bx, int cx)
    {
        CellKey const key{cx};
        if (auto it = cache.find(key); it != cache.end())
            return it->second;
        double const d = bx(cx + 1) - bx(cx);
        cache.emplace(key, d);
        return d;
    }

    static double subzoneDiv2d_(auto& cache, auto& bx, auto& by, int cx, int cy)
    {
        CellKey const key{cx, cy};
        if (auto it = cache.find(key); it != cache.end())
            return it->second;
        double const d = (bx(cx + 1, cy) - bx(cx, cy)) + (by(cx, cy + 1) - by(cx, cy));
        cache.emplace(key, d);
        return d;
    }

    static double subzoneDiv3d_(auto& cache, auto& bx, auto& by, auto& bz, int cx, int cy, int cz)
    {
        CellKey const key{cx, cy, cz};
        if (auto it = cache.find(key); it != cache.end())
            return it->second;
        double const d = (bx(cx + 1, cy, cz) - bx(cx, cy, cz))
                         + (by(cx, cy + 1, cz) - by(cx, cy, cz))
                         + (bz(cx, cy, cz + 1) - bz(cx, cy, cz));
        cache.emplace(key, d);
        return d;
    }

    ResMan& rm_;
    int b_id_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_ADPT_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
