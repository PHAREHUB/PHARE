#ifndef PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "amr/resources_manager/amr_utils.hpp"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "core/utilities/constants.hpp"

#include <cassert>
#include <iostream>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

template<typename ResMan, typename FieldDataT>
class MagneticRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    using Geometry        = typename FieldDataT::Geometry;
    using gridlayout_type = typename FieldDataT::gridlayout_type;

    static constexpr std::size_t dimension = FieldDataT::dimension;

    MagneticRefinePatchStrategy(ResMan& resourcesManager,
                                std::optional<std::shared_ptr<SAMRAI::hier::PatchLevel>> old_level
                                = std::nullopt)
        : rm_{resourcesManager}
        , oldLevel_{std::move(old_level)}
    {
    }

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                       const SAMRAI::hier::IntVector& ghost_width_to_fill)
    {
    }

    SAMRAI::hier::IntVector getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const
    {
        return SAMRAI::hier::IntVector(dim, 1); // hard-coded 0th order base interpolation
    }


    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box, SAMRAI::hier::IntVector const& ratio)
    {
    }

    // We compute the values of the new fine magnetic faces using what was already refined, ie the
    // values on the old coarse faces.
    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box, SAMRAI::hier::IntVector const& ratio)
    {
        auto bx_id = rm_.getID("EM_B_x");
        auto by_id = rm_.getID("EM_B_y");
        auto bz_id = rm_.getID("EM_B_z");

        auto& bx = FieldDataT::getField(fine, *bx_id);
        auto& by = FieldDataT::getField(fine, *by_id);
        auto& bz = FieldDataT::getField(fine, *bz_id);

        auto layout        = PHARE::amr::layoutFromPatch<gridlayout_type>(fine);
        auto fineBoxLayout = Geometry::layoutFromBox(fine_box, layout);

        SAMRAI::hier::Box fine_box_x
            = Geometry::toFieldBox(fine_box, bx.physicalQuantity(), fineBoxLayout, false);
        SAMRAI::hier::Box fine_box_y
            = Geometry::toFieldBox(fine_box, by.physicalQuantity(), fineBoxLayout, false);
        SAMRAI::hier::Box fine_box_z
            = Geometry::toFieldBox(fine_box, bz.physicalQuantity(), fineBoxLayout, false);

        if (oldLevel_)
        {
            SAMRAI::hier::BoxContainer oldBoxes = (*oldLevel_)->getBoxes();

            SAMRAI::hier::BoxContainer oldCoarseFineIntersectionX
                = Geometry::toFieldBoxes(oldBoxes, bx.physicalQuantity(), layout, false);
            oldCoarseFineIntersectionX.intersectBoxes(fine_box_x);

            SAMRAI::hier::BoxContainer oldCoarseFineIntersectionY
                = Geometry::toFieldBoxes(oldBoxes, by.physicalQuantity(), layout, false);
            oldCoarseFineIntersectionY.intersectBoxes(fine_box_y);

            SAMRAI::hier::BoxContainer oldCoarseFineIntersectionZ
                = Geometry::toFieldBoxes(oldBoxes, bz.physicalQuantity(), layout, false);
            oldCoarseFineIntersectionZ.intersectBoxes(fine_box_z);

            auto imposeCopiedOldCoarseFine = [&](auto& b, auto const& id,
                                                 auto const& oldCFintersection) {
                for (auto const& box : oldCFintersection)
                {
                    auto oldPatch  = (*oldLevel_)->getPatch(box.getBoxId());
                    auto& bOld     = FieldDataT::getField(*oldPatch, id);
                    auto layoutOld = PHARE::amr::layoutFromPatch<gridlayout_type>(*oldPatch);

                    // find the old box so that we can get the local index for it
                    // maybe doable with just the layout?
                    // SAMRAI::hier::BoxContainer old_overlap_boxes;
                    // oldBoxes.findOverlapBoxes(old_overlap_boxes, box);
                    // assert(old_overlap_boxes.size() == 1);
                    // SAMRAI::hier::Box old_box = old_overlap_boxes.front();

                    int iStartX = box.lower(dirX);
                    int iEndX   = box.upper(dirX);

                    if constexpr (dimension == 1)
                    {
                        for (int amrix = iStartX; amrix <= iEndX; ++amrix)
                        {
                            auto idxf   = layout.AMRToLocal(core::Point{amrix});
                            auto idxold = layoutOld.AMRToLocal(core::Point{amrix});
                            auto ixf    = idxf[dirX];
                            auto ixold  = idxold[dirX];
                            b(ixf)      = bOld(ixold);
                        }
                    }
                    else if constexpr (dimension >= 2)
                    {
                        int iStartY = box.lower(dirY);
                        int iEndY   = box.upper(dirY);

                        if constexpr (dimension == 2)
                        {
                            for (int amrix = iStartX; amrix <= iEndX; ++amrix)
                            {
                                for (int amriy = iStartY; amriy <= iEndY; ++amriy)
                                {
                                    auto idxf   = layout.AMRToLocal(core::Point{amrix, amriy});
                                    auto idxold = layoutOld.AMRToLocal(core::Point{amrix, amriy});
                                    auto ixf    = idxf[dirX];
                                    auto iyf    = idxf[dirY];
                                    auto ixold  = idxold[dirX];
                                    auto iyold  = idxold[dirY];
                                    b(ixf, iyf) = bOld(ixold, iyold);
                                }
                            }
                        }
                        else if constexpr (dimension == 3)
                        {
                            int iStartZ = box.lower(dirZ);
                            int iEndZ   = box.upper(dirZ);

                            for (int amriz = iStartZ; amriz <= iEndZ; ++amriz)
                            {
                                for (int amriy = iStartY; amriy <= iEndY; ++amriy)
                                {
                                    for (int amrix = iStartX; amrix <= iEndX; ++amrix)
                                    {
                                        auto idxf
                                            = layout.AMRToLocal(core::Point{amrix, amriy, amriz});
                                        auto idxold = layoutOld.AMRToLocal(
                                            core::Point{amrix, amriy, amriz});
                                        auto ixf         = idxf[dirX];
                                        auto iyf         = idxf[dirY];
                                        auto izf         = idxf[dirZ];
                                        auto ixold       = idxold[dirX];
                                        auto iyold       = idxold[dirY];
                                        auto izold       = idxold[dirZ];
                                        b(ixf, iyf, izf) = bOld(ixold, iyold, izold);
                                    }
                                }
                            }
                        }
                    }
                }
            };

            imposeCopiedOldCoarseFine(bx, *bx_id, oldCoarseFineIntersectionX);
            imposeCopiedOldCoarseFine(by, *by_id, oldCoarseFineIntersectionY);
            imposeCopiedOldCoarseFine(bz, *bz_id, oldCoarseFineIntersectionZ);
        }

        int ixStartX = fine_box_x.lower(dirX);
        int ixEndX   = fine_box_x.upper(dirX);

        if constexpr (dimension == 1)
        {
            for (int amrix = ixStartX; amrix <= ixEndX; ++amrix)
            {
                auto i  = layout.AMRToLocal(core::Point{amrix});
                auto ix = i[dirX];

                if (ix % 2 == 1)
                    bx(ix) = 0.5 * (bx(ix - 1) + bx(ix + 1));
            }
        }

        else if constexpr (dimension >= 2)
        {
            int iyStartX = fine_box_x.lower(dirY);
            int iyEndX   = fine_box_x.upper(dirY);
            int ixStartY = fine_box_y.lower(dirX);
            int ixEndY   = fine_box_y.upper(dirX);
            int iyStartY = fine_box_y.lower(dirY);
            int iyEndY   = fine_box_y.upper(dirY);

            auto p_plus  = [](auto const i, auto const offset) { return i + 2 - offset; };
            auto p_minus = [](auto const i, auto const offset) { return i - offset; };

            auto d_plus  = [](auto const i, auto const offset) { return i + 1 - offset; };
            auto d_minus = [](auto const i, auto const offset) { return i - offset; };

            if constexpr (dimension == 2)
            {
                for (int amrix = ixStartX; amrix <= ixEndX; ++amrix)
                {
                    for (int amriy = iyStartX; amriy <= iyEndX; ++amriy)
                    {
                        auto i  = layout.AMRToLocal(core::Point{amrix, amriy});
                        auto ix = i[dirX];
                        auto iy = i[dirY];
                        //                            | <- here with offset = 1
                        //                          -- --
                        //                            | <- or here with offset = 0
                        if (ix % 2 == 1)
                        {
                            // If dual no offset, ie primal for the field we are actually
                            // modifying, but dual for the field we are indexing to compute
                            // second and third order terms, then the formula reduces to offset
                            // = 1
                            int xoffset = 1;
                            int yoffset = (iy % 2 == 0) ? 0 : 1;

                            bx(ix, iy) = 0.5 * (bx(ix - 1, iy) + bx(ix + 1, iy))
                                         + 0.25
                                               * (by(d_minus(ix, xoffset), p_minus(iy, yoffset))
                                                  - by(d_minus(ix, xoffset), p_plus(iy, yoffset))
                                                  - by(d_plus(ix, xoffset), p_minus(iy, yoffset))
                                                  + by(d_plus(ix, xoffset), p_plus(iy, yoffset)));
                        }
                    }
                }

                for (int amrix = ixStartY; amrix <= ixEndY; ++amrix)
                {
                    for (int amriy = iyStartY; amriy <= iyEndY; ++amriy)
                    {
                        auto i  = layout.AMRToLocal(core::Point{amrix, amriy});
                        auto ix = i[dirX];
                        auto iy = i[dirY];
                        //                            |
                        //  here with offset = 0 -> -- -- <- or here with offset = 1
                        //                            |
                        if (iy % 2 == 1)
                        {
                            int xoffset = (ix % 2 == 0) ? 0 : 1;
                            int yoffset = 1;

                            by(ix, iy) = 0.5 * (by(ix, iy - 1) + by(ix, iy + 1))
                                         + 0.25
                                               * (bx(p_minus(ix, xoffset), d_minus(iy, yoffset))
                                                  - bx(p_plus(ix, xoffset), d_minus(iy, yoffset))
                                                  - bx(p_minus(ix, xoffset), d_plus(iy, yoffset))
                                                  + bx(p_plus(ix, xoffset), d_plus(iy, yoffset)));
                        }
                    }
                }
            }

            else if constexpr (dimension == 3)
            {
                int izStartX = fine_box_x.lower(dirZ);
                int izEndX   = fine_box_x.upper(dirZ);
                int izStartY = fine_box_y.lower(dirZ);
                int izEndY   = fine_box_y.upper(dirZ);
                int ixStartZ = fine_box_z.lower(dirX);
                int ixEndZ   = fine_box_z.upper(dirX);
                int iyStartZ = fine_box_z.lower(dirY);
                int iyEndZ   = fine_box_z.upper(dirY);
                int izStartZ = fine_box_z.lower(dirZ);
                int izEndZ   = fine_box_z.upper(dirZ);

                auto Dx = layout.meshSize()[dirX];
                auto Dy = layout.meshSize()[dirY];
                auto Dz = layout.meshSize()[dirZ];

                // Toth and Roe (2002) use a formulation we the indexing is centered
                // on the coarse cell. Since this is not our case, we need to have a
                // different offset for indexing and applying the +-1 factor to the
                // third order terms. That's the job of the ijk_factor lambda
                auto ijk_factor = [](auto const offset) { return offset == 0 ? -1 : 1; };

                for (int amrix = ixStartX; amrix <= ixEndX; ++amrix)
                {
                    for (int amriy = iyStartX; amriy <= iyEndX; ++amriy)
                    {
                        for (int amriz = izStartX; amriz <= izEndX; ++amriz)
                        {
                            auto i  = layout.AMRToLocal(core::Point{amrix, amriy, amriz});
                            auto ix = i[dirX];
                            auto iy = i[dirY];
                            auto iz = i[dirZ];


                            if (ix % 2 == 1)
                            {
                                int xoffset = 1;
                                int yoffset = (iy % 2 == 0) ? 0 : 1;
                                int zoffset = (iz % 2 == 0) ? 0 : 1;

                                bx(ix, iy, iz)
                                    = 0.5 * (bx(ix - 1, iy, iz) + bx(ix + 1, iy, iz))
                                      + 0.125
                                            * (by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                  d_minus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset)))
                                      + 0.125
                                            * (bz(d_minus(ix, xoffset), d_minus(iy, yoffset),
                                                  p_minus(iz, zoffset))
                                               + bz(d_minus(ix, xoffset), d_plus(iy, yoffset),
                                                    p_minus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), d_minus(iy, yoffset),
                                                    p_minus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), d_plus(iy, yoffset),
                                                    p_minus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), d_minus(iy, yoffset),
                                                    p_plus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), d_plus(iy, yoffset),
                                                    p_plus(iz, zoffset))
                                               + bz(d_plus(ix, xoffset), d_minus(iy, yoffset),
                                                    p_plus(iz, zoffset))
                                               + bz(d_plus(ix, xoffset), d_plus(iy, yoffset),
                                                    p_plus(iz, zoffset)))
                                      + (0.125 * ijk_factor(zoffset) * Dz * Dz
                                         / (Dx * Dx + Dz * Dz))
                                            * (by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                  d_plus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset)))
                                      + (0.125 * ijk_factor(yoffset) * Dy * Dy
                                         / (Dx * Dx + Dy * Dy))
                                            * (bz(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                  d_plus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset)));
                            }
                        }
                    }
                }

                for (int amrix = ixStartY; amrix <= ixEndY; ++amrix)
                {
                    for (int amriy = iyStartY; amriy <= iyEndY; ++amriy)
                    {
                        for (int amriz = izStartY; amriz <= izEndY; ++amriz)
                        {
                            auto i  = layout.AMRToLocal(core::Point{amrix, amriy, amriz});
                            auto ix = i[dirX];
                            auto iy = i[dirY];
                            auto iz = i[dirZ];

                            if (iy % 2 == 1)
                            {
                                int xoffset = (ix % 2 == 0) ? 0 : 1;
                                int yoffset = 1;
                                int zoffset = (iz % 2 == 0) ? 0 : 1;

                                by(ix, iy, iz)
                                    = 0.5 * (by(ix, iy - 1, iz) + by(ix, iy + 1, iz))
                                      + 0.125
                                            * (bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                  d_minus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset)))
                                      + 0.125
                                            * (bz(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                  d_minus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + bz(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + bz(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset)))
                                      + (0.125 * ijk_factor(xoffset) * Dx * Dx
                                         / (Dx * Dx + Dy * Dy))
                                            * (bz(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                  d_plus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bz(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bz(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset)))
                                      + (0.125 * ijk_factor(zoffset) * Dz * Dz
                                         / (Dy * Dy + Dz * Dz))
                                            * (bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                  d_plus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset)));
                            }
                        }
                    }
                }

                for (int amrix = ixStartZ; amrix <= ixEndZ; ++amrix)
                {
                    for (int amriy = iyStartZ; amriy <= iyEndZ; ++amriy)
                    {
                        for (int amriz = izStartZ; amriz <= izEndZ; ++amriz)
                        {
                            auto i  = layout.AMRToLocal(core::Point{amrix, amriy, amriz});
                            auto ix = i[dirX];
                            auto iy = i[dirY];
                            auto iz = i[dirZ];
                            if (iz % 2 == 1)
                            {
                                int xoffset = (ix % 2 == 0) ? 0 : 1;
                                int yoffset = (iy % 2 == 0) ? 0 : 1;
                                int zoffset = 1;

                                bz(ix, iy, iz)
                                    = 0.5 * (bz(ix, iy, iz - 1) + bz(ix, iy, iz + 1))
                                      + 0.125
                                            * (bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                  d_minus(iz, zoffset))
                                               + bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset)))
                                      + 0.125
                                            * (by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                  d_minus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               + by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset)))
                                      + (0.125 * ijk_factor(yoffset) * Dy * Dy
                                         / (Dy * Dy + Dz * Dz))
                                            * (bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                  d_plus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - bx(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset)))
                                      + (0.125 * ijk_factor(xoffset) * Dx * Dx
                                         / (Dx * Dx + Dz * Dz))
                                            * (by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                  d_plus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_plus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_plus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_minus(ix, xoffset), p_plus(iy, yoffset),
                                                    d_minus(iz, zoffset))
                                               + by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_plus(iz, zoffset))
                                               - by(d_minus(ix, xoffset), p_minus(iy, yoffset),
                                                    d_minus(iz, zoffset)));
                            }
                        }
                    }
                }
            }
        }
    }

private:
    ResMan& rm_;
    std::optional<std::shared_ptr<SAMRAI::hier::PatchLevel>> oldLevel_;
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
