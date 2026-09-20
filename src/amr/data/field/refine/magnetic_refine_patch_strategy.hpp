#ifndef PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/constants.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/refine/coarse_cell_round_out.hpp"

#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <array>
#include <cassert>
#include <optional>
#include <stdexcept>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

// PHARE::core::Box counterpart of coarse_cell_round_out.hpp's SAMRAI::hier::Box helpers - same
// parity rules (roundDownToEven/roundUpToOddIndex/roundUpToEvenIndex/roundDownToOddIndex/
// isOddIndex), applied to the per-component field boxes tiles are natively expressed in, to avoid
// round-tripping through SAMRAI::hier::Box per tile.
//
// Shrinks box to the largest whole-coarse-cell union it contains, per direction - used to clip a
// reconstruction region down to what a single tile's own ghost_box() can safely process, instead
// of rejecting the tile outright when the region only partially fits (a tile with any overlap
// should still handle the safe part of it, since other overlapping tiles hold their own separate
// copies of the shared cells and each needs its own copy written).
template<std::size_t dim>
NO_DISCARD std::optional<core::Box<int, dim>>
roundFieldBoxInToCoarseCells(core::Box<int, dim> box,
                             std::array<core::QtyCentering, dim> const& centering)
{
    for (std::size_t d = 0; d < dim; ++d)
    {
        box.lower[d] = roundUpToEvenIndex(box.lower[d]);
        box.upper[d] = centering[d] == core::QtyCentering::primal ? roundDownToEven(box.upper[d])
                                                                  : roundDownToOddIndex(box.upper[d]);
        if (box.lower[d] > box.upper[d])
            return std::nullopt;
    }
    return box;
}


template<typename ResMan, typename TensorFieldDataT>
class MagneticRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
    auto make_fine_field_boxes(auto& fields, auto& fine_box, auto& layout, auto& fineLayout) const
    {
        return core::for_N_make_array<N>([&](auto i) {
            using PhysicalQuantity = std::decay_t<decltype(fields[i].physicalQuantity())>;

            return phare_box_from<dimension>(
                FieldGeometry<gridlayout_type, PhysicalQuantity>::toFieldBox(
                    fine_box, fields[i].physicalQuantity(), fineLayout));
        });
    }

public:
    using Geometry        = TensorFieldDataT::Geometry;
    using gridlayout_type = TensorFieldDataT::gridlayout_type;

    static constexpr std::size_t N         = TensorFieldDataT::N;
    static constexpr std::size_t dimension = TensorFieldDataT::dimension;

    MagneticRefinePatchStrategy(ResMan& resourcesManager)
        : rm_{resourcesManager}
        , b_id_{-1}
    {
    }

    void assertIDsSet() const
    {
        assert(b_id_ >= 0 && "MagneticRefinePatchStrategy: IDs must be registered before use");
    }

    void registerIDs(int const b_id) { b_id_ = b_id; }

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                       SAMRAI::hier::IntVector const& ghost_width_to_fill) override
    {
    }

    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector(dim, 1); // hard-coded 0th order base interpolation
    }


    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box,
                          SAMRAI::hier::IntVector const& ratio) override
    {
    }

    // The region a magnetic reconstruction over the SAMRAI fill box `fine_box` must run over.
    //
    // The Toth & Roe reconstruction below (fix()/postprocessBx3d & friends) touches faces of
    // exactly one coarse cell per "new fine face" it fills. SAMRAI's fill boxes do not respect
    // that: a recursive refine schedule fills a coarse-interpolation temporary plus a 1-cell
    // stencil ring, and one cell is always *half* a coarse cell, whatever the temporary's own
    // parity - so the coarse-aligned face on the far side of a cut cell may never have been
    // refined yet. Rounding the box out to a whole-coarse-cell union (see
    // coarse_cell_round_out.hpp) guarantees every input the reconstruction reads was written.
    // See https://github.com/PHAREHUB/PHARE/issues/1296.
    SAMRAI::hier::Box reconstructionRegion(SAMRAI::hier::Box const& fine_box,
                                           SAMRAI::hier::Box const& ghost_box) const
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

    // Tile-level analog of reconstructionRegion, in already-per-component field-box space.
    //
    // fine_field_box is already a whole-coarse-cell union at the PATCH level, but tiles are not
    // always coarse-cell aligned (refined-level patches can have odd extents, which pushes their
    // tiling off the coarse grid - see tile_set_mapper.hpp's all_even_shape fallback), so clipping
    // it to one tile's ghost_box() can still cut a coarse cell in half at the tile boundary.
    //
    // Field ghost cells have no single owning tile (unlike particles): tiles' ghost boxes overlap,
    // and each tile holds its own separate copy of any cell in that overlap, so every tile with
    // some overlap needs its own copy written, not just whichever one is asked first. Rounding the
    // clipped overlap IN (shrinking to the largest whole-coarse-cell union that fits) rather than
    // OUT (growing, then rejecting the whole thing if it no longer fits back in the tile) means a
    // tile only ever gets asked to do what it can safely do with its own stored data, and still
    // contributes that safe subset instead of being skipped entirely over a sliver at its edge.
    static std::optional<core::Box<int, dimension>>
    tileReconstructionFieldBox(core::Box<int, dimension> const& fine_field_box,
                               core::Box<int, dimension> const& tile_ghost_box,
                               std::array<core::QtyCentering, dimension> const& centering)
    {
        auto const clipped = fine_field_box * tile_ghost_box;
        if (!clipped)
            return std::nullopt;

        return roundFieldBoxInToCoarseCells<dimension>(*clipped, centering);
    }

    // We compute the values of the new fine magnetic faces using what was already refined, ie
    // the values on the old coarse faces.
    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
        assertIDsSet();

        auto& fields       = TensorFieldDataT::getFields(fine, b_id_);
        auto& [bx, by, bz] = fields;

        auto const region
            = reconstructionRegion(fine_box, fine.getPatchData(b_id_)->getGhostBox());

        auto layout                 = PHARE::amr::layoutFromPatch<gridlayout_type>(fine);
        auto fineBoxLayout          = Geometry::layoutFromBox(region, layout);
        auto const fine_field_boxes = make_fine_field_boxes(fields, region, layout, fineBoxLayout);

        using Field_t = std::decay_t<decltype(bx)>;
        if constexpr (core::is_field_tile_set_v<Field_t>)
        {
            auto const centerings = core::for_N_make_array<N>(
                [&](auto i) { return gridlayout_type::centering(fields[i].physicalQuantity()); });

            for (std::size_t ti = 0; ti < bx().size(); ++ti)
            {
                auto& bx_tile = bx()[ti];
                auto& by_tile = by()[ti];
                auto& bz_tile = bz()[ti];

                auto const xoverlap = tileReconstructionFieldBox(
                    fine_field_boxes[dirX], bx_tile.ghost_box(), centerings[dirX]);
                auto const yoverlap = tileReconstructionFieldBox(
                    fine_field_boxes[dirY], by_tile.ghost_box(), centerings[dirY]);
                auto const zoverlap = tileReconstructionFieldBox(
                    fine_field_boxes[dirZ], bz_tile.ghost_box(), centerings[dirZ]);

                if (!xoverlap && !yoverlap && !zoverlap)
                    continue;

                auto const& tile_layout = bx_tile.layout();
                fix(bx_tile(), by_tile(), bz_tile(), tile_layout,
                    std::array{xoverlap, yoverlap, zoverlap});
            }
        }
        else
        {
            fix(bx, by, bz, layout,
                std::array<std::optional<core::Box<int, dimension>>, N>{
                    fine_field_boxes[dirX], fine_field_boxes[dirY], fine_field_boxes[dirZ]});
        }
    }


    void fix(auto& bx, auto& by, auto& bz, auto& layout, auto const& fine_field_boxes)
    {
        if constexpr (dimension == 1)
        {
            // if we ever go to c++23 we could use std::views::zip to iterate both on the local and
            // global indices instead of passing the box to do an amr to local inside the function,
            // which is not obvious at call site
            if (fine_field_boxes[dirX])
                for (auto const& i : *fine_field_boxes[dirX])
                    postprocessBx1d(bx, layout, i);
        }

        else if constexpr (dimension == 2)
        {
            if (fine_field_boxes[dirX])
                for (auto const& i : *fine_field_boxes[dirX])
                    postprocessBx2d(bx, by, layout, i);


            if (fine_field_boxes[dirY])
                for (auto const& i : *fine_field_boxes[dirY])
                    postprocessBy2d(bx, by, layout, i);
        }

        else if constexpr (dimension == 3)
        {
            auto meshSize = layout.meshSize();

            if (fine_field_boxes[dirX])
                for (auto const& i : *fine_field_boxes[dirX])
                    postprocessBx3d(bx, by, bz, meshSize, layout, i);


            if (fine_field_boxes[dirY])
                for (auto const& i : *fine_field_boxes[dirY])
                    postprocessBy3d(bx, by, bz, meshSize, layout, i);


            if (fine_field_boxes[dirZ])
                for (auto const& i : *fine_field_boxes[dirZ])
                    postprocessBz3d(bx, by, bz, meshSize, layout, i);
        }
    }




    static auto isNewFineFace(auto const& amrIdx, auto const dir)
    {
        // amr index can be negative so test !=0 and not ==1
        // to see if this is odd or even
        return amrIdx[dir] % 2 != 0;
    }

    static void postprocessBx1d(auto& bx, auto const& layout, core::Point<int, dimension> idx)
    {
        auto const locIdx = layout.AMRToLocal(idx);
        auto const ix     = locIdx[dirX];
        if (isNewFineFace(idx, dirX))
            bx(ix) = 0.5 * (bx(ix - 1) + bx(ix + 1));
    }

    static void postprocessBx2d(auto& bx, auto& by, auto const& layout,
                                core::Point<int, dimension> idx)
    {
        auto const locIdx = layout.AMRToLocal(idx);
        auto const ix     = locIdx[dirX];
        auto const iy     = locIdx[dirY];
        //                            | <- here with offset = 1
        //                          -- --
        //                            | <- or here with offset = 0
        if (isNewFineFace(idx, dirX))
        {
            // If dual no offset, ie primal for the field we are actually
            // modifying, but dual for the field we are indexing to compute
            // second and third order terms, then the formula reduces to offset
            // = 1
            int const xoffset = 1;
            int const yoffset = (idx[dirY] % 2 == 0) ? 0 : 1;

            bx(ix, iy) = 0.5 * (bx(ix - 1, iy) + bx(ix + 1, iy))
                         + 0.25
                               * (by(d_minus(ix, xoffset), p_minus(iy, yoffset))
                                  - by(d_minus(ix, xoffset), p_plus(iy, yoffset))
                                  - by(d_plus(ix, xoffset), p_minus(iy, yoffset))
                                  + by(d_plus(ix, xoffset), p_plus(iy, yoffset)));
        }
    }

    static void postprocessBy2d(auto& bx, auto& by, auto const& layout,
                                core::Point<int, dimension> idx)
    {
        auto const locIdx = layout.AMRToLocal(idx);
        auto const ix     = locIdx[dirX];
        auto const iy     = locIdx[dirY];
        //                            |
        //  here with offset = 0 -> -- -- <- or here with offset = 1
        //                            |
        if (isNewFineFace(idx, dirY))
        {
            int const xoffset = (idx[dirX] % 2 == 0) ? 0 : 1;
            int const yoffset = 1;

            by(ix, iy) = 0.5 * (by(ix, iy - 1) + by(ix, iy + 1))
                         + 0.25
                               * (bx(p_minus(ix, xoffset), d_minus(iy, yoffset))
                                  - bx(p_plus(ix, xoffset), d_minus(iy, yoffset))
                                  - bx(p_minus(ix, xoffset), d_plus(iy, yoffset))
                                  + bx(p_plus(ix, xoffset), d_plus(iy, yoffset)));
        }
    }

    static void postprocessBx3d(auto& bx, auto& by, auto& bz, auto const& meshSize,
                                auto const& layout, core::Point<int, dimension> idx)
    {
        auto const Dx = meshSize[dirX];
        auto const Dy = meshSize[dirY];
        auto const Dz = meshSize[dirZ];

        auto const locIdx = layout.AMRToLocal(idx);
        auto const ix     = locIdx[dirX];
        auto const iy     = locIdx[dirY];
        auto const iz     = locIdx[dirZ];

        if (isNewFineFace(idx, dirX))
        {
            int const xoffset = 1;
            int const yoffset = (idx[dirY] % 2 == 0) ? 0 : 1;
            int const zoffset = (idx[dirZ] % 2 == 0) ? 0 : 1;

            bx(ix, iy, iz)
                = 0.5 * (bx(ix - 1, iy, iz) + bx(ix + 1, iy, iz))
                  + 0.125
                        * (by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           + by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset)))
                  + 0.125
                        * (bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset))
                           + bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           + bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset)))
                  + (0.125 * ijk_factor_[zoffset] * Dz * Dz / (Dx * Dx + Dz * Dz))
                        * (by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset)))
                  + (0.125 * ijk_factor_[yoffset] * Dy * Dy / (Dx * Dx + Dy * Dy))
                        * (bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset)));
        }
    };

    static void postprocessBy3d(auto& bx, auto& by, auto& bz, auto const& meshSize,
                                auto const& layout, core::Point<int, dimension> idx)
    {
        auto const Dx = meshSize[dirX];
        auto const Dy = meshSize[dirY];
        auto const Dz = meshSize[dirZ];

        auto const locIdx = layout.AMRToLocal(idx);
        auto const ix     = locIdx[dirX];
        auto const iy     = locIdx[dirY];
        auto const iz     = locIdx[dirZ];

        if (isNewFineFace(idx, dirY))
        {
            int const xoffset = (idx[dirX] % 2 == 0) ? 0 : 1;
            int const yoffset = 1;
            int const zoffset = (idx[dirZ] % 2 == 0) ? 0 : 1;

            by(ix, iy, iz)
                = 0.5 * (by(ix, iy - 1, iz) + by(ix, iy + 1, iz))
                  + 0.125
                        * (bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           + bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset)))
                  + 0.125
                        * (bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           + bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           + bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset)))
                  + (0.125 * ijk_factor_[xoffset] * Dx * Dx / (Dx * Dx + Dy * Dy))
                        * (bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_plus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_plus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_minus(ix, xoffset), d_plus(iy, yoffset), p_minus(iz, zoffset))
                           + bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_plus(iz, zoffset))
                           - bz(d_minus(ix, xoffset), d_minus(iy, yoffset), p_minus(iz, zoffset)))
                  + (0.125 * ijk_factor_[zoffset] * Dz * Dz / (Dy * Dy + Dz * Dz))
                        * (bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset)));
        }
    };

    static void postprocessBz3d(auto& bx, auto& by, auto& bz, auto const& meshSize,
                                auto const& layout, core::Point<int, dimension> idx)
    {
        auto const Dx = meshSize[dirX];
        auto const Dy = meshSize[dirY];
        auto const Dz = meshSize[dirZ];

        auto const locIdx = layout.AMRToLocal(idx);
        auto const ix     = locIdx[dirX];
        auto const iy     = locIdx[dirY];
        auto const iz     = locIdx[dirZ];

        if (isNewFineFace(idx, dirZ))
        {
            int const xoffset = (idx[dirX] % 2 == 0) ? 0 : 1;
            int const yoffset = (idx[dirY] % 2 == 0) ? 0 : 1;
            int const zoffset = 1;

            bz(ix, iy, iz)
                = 0.5 * (bz(ix, iy, iz - 1) + bz(ix, iy, iz + 1))
                  + 0.125
                        * (bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset))
                           + bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           + bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset)))
                  + 0.125
                        * (by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           + by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           + by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset)))
                  + (0.125 * ijk_factor_[yoffset] * Dy * Dy / (Dy * Dy + Dz * Dz))
                        * (bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_plus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_plus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_minus(ix, xoffset), d_plus(iy, yoffset), d_minus(iz, zoffset))
                           + bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_plus(iz, zoffset))
                           - bx(p_minus(ix, xoffset), d_minus(iy, yoffset), d_minus(iz, zoffset)))
                  + (0.125 * ijk_factor_[xoffset] * Dx * Dx / (Dx * Dx + Dz * Dz))
                        * (by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_plus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_plus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_minus(ix, xoffset), p_plus(iy, yoffset), d_minus(iz, zoffset))
                           + by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_plus(iz, zoffset))
                           - by(d_minus(ix, xoffset), p_minus(iy, yoffset), d_minus(iz, zoffset)));
        }
    };



private:
    static auto p_plus(auto const i, auto const offset) { return i + 2 - offset; };
    static auto p_minus(auto const i, auto const offset) { return i - offset; };

    static auto d_plus(auto const i, auto const offset) { return i + 1 - offset; };
    static auto d_minus(auto const i, auto const offset) { return i - offset; };

    // Toth and Roe (2002) use a formulation we the indexing is centered
    // on the coarse cell. Since this is not our case, we need to have a
    // different offset for indexing and applying the +-1 factor to the
    // third order terms. That's the job of the ijk_factor_ array.
    static constexpr std::array<int, 2> ijk_factor_{-1, 1};

    ResMan& rm_;
    int b_id_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
