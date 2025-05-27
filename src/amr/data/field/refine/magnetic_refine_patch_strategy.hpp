#ifndef PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/utilities/constants.hpp"

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

template<typename ResMan, typename FieldDataT>
class MagneticRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    MagneticRefinePatchStrategy(ResMan& resourcesManager)
        : rm_{resourcesManager}
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

        int iStartX = fine_box.lower(dirX);
        int iEndX   = fine_box.upper(dirX);

        if constexpr (FieldDataT::dimension == 1)
        {
            for (int ix = iStartX; ix <= iEndX; ++ix)
            {
                if (ix % 2 == 1)
                    bx(ix) = 0.5 * (bx(ix - 1) + bx(ix + 1));
            }
        }

        else if constexpr (FieldDataT::dimension >= 2)
        {
            int iStartY = fine_box.lower(dirY);
            int iEndY   = fine_box.upper(dirY);

            auto p_plus  = [](auto const i, auto const offset) { return i + 2 - offset; };
            auto p_minus = [](auto const i, auto const offset) { return i - offset; };

            auto d_plus  = [](auto const i, auto const offset) { return i + 1 - offset; };
            auto d_minus = [](auto const i, auto const offset) { return i - offset; };

            if constexpr (FieldDataT::dimension == 2)
            {
                for (int ix = iStartX; ix <= iEndX; ++ix)
                {
                    for (int iy = iStartY; iy <= iEndY; ++iy)
                    {
                        //                            | <- here with offset = 1
                        //                          -- --
                        //                            | <- or here with offset = 0
                        if (ix % 2 == 1)
                        {
                            // If dual no offset, ie primal for the field we are actually modifying,
                            // but dual for the field we are indexing to compute second and third
                            // order terms, then the formula reduces to offset = 1
                            int xoffset = 1;
                            int yoffset = (iy % 2 == 0) ? 0 : 1;

                            bx(ix, iy) = 0.5 * (bx(ix - 1, iy) + bx(ix + 1, iy))
                                         + 0.25
                                               * (by(d_minus(ix, xoffset), p_minus(iy, yoffset))
                                                  - by(d_minus(ix, xoffset), p_plus(iy, yoffset))
                                                  - by(d_plus(ix, xoffset), p_minus(iy, yoffset))
                                                  + by(d_plus(ix, xoffset), p_plus(iy, yoffset)));
                        }
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

            else if constexpr (FieldDataT::dimension == 3)
            {
                int iStartZ = fine_box.lower(dirZ);
                int iEndZ   = fine_box.upper(dirZ);

                auto layout = PHARE::amr::layoutFromPatch<FieldDataT::gridlayout_type>(fine);

                auto Dx = layout.meshSize()[dirX];
                auto Dy = layout.meshSize()[dirY];
                auto Dz = layout.meshSize()[dirZ];

                for (int ix = iStartX; ix <= iEndX; ++ix)
                {
                    for (int iy = iStartY; iy <= iEndY; ++iy)
                    {
                        for (int iz = iStartZ; iz <= iEndZ; ++iz)
                        {
                            // Toth and Roe (2002) use a formulation we the indexing is centered
                            // on the coarse cell. Since this is not our case, we need to have a
                            // different offset for indexing and applying the +-1 factor to the
                            // third order terms. That's the job of the ijk_factor lambda
                            auto ijk_factor
                                = [](auto const offset) { return offset == 0 ? -1 : 1; };

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

                            if (iz % 2 == 1)
                            {
                                int xoffset = (ix % 2 == 0) ? 0 : 1;
                                int yoffset = (iy % 2 == 0) ? 0 : 1;
                                int zoffset = 1;

                                bz(ix, iy, iz)
                                    = 0.5 * (bz(ix, iy, iz - 1) + bz(ix, iy, iz - 1))
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
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
