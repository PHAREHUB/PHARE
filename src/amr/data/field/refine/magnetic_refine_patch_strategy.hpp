#ifndef PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "amr/data/field/field_geometry.hpp"
#include "core/utilities/constants.hpp"


#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/resources_manager/amr_utils.hpp"


#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "core/utilities/types.hpp"

#include <array>
#include <cassert>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

template<typename ResMan, typename TensorFieldDataT>
class MagneticRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    using Geometry        = typename TensorFieldDataT::Geometry;
    using gridlayout_type = typename TensorFieldDataT::gridlayout_type;

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
                                       const SAMRAI::hier::IntVector& ghost_width_to_fill) override
    {
    }

    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(const SAMRAI::tbox::Dimension& dim) const override
    {
        return SAMRAI::hier::IntVector(dim, 1); // hard-coded 0th order base interpolation
    }


    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box,
                          SAMRAI::hier::IntVector const& ratio) override
    {
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

        auto layout        = PHARE::amr::layoutFromPatch<gridlayout_type>(fine);
        auto fineBoxLayout = Geometry::layoutFromBox(fine_box, layout);

        auto fine_field_box = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            using PhysicalQuantity = std::decay_t<decltype(fields[i].physicalQuantity())>;

            return FieldGeometry<gridlayout_type, PhysicalQuantity>::toFieldBox(
                fine_box, fields[i].physicalQuantity(), fineBoxLayout);
        });

        if constexpr (dimension == 1)
        {
            // if we ever go to c++23 we could use std::views::zip to iterate both on the local and
            // global indices instead of passing the box to do an amr to local inside the function,
            // which is not obvious at call site
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
            {
                postprocessBx1d(bx, layout, i);
            }
        }

        else if constexpr (dimension == 2)
        {
            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
            {
                postprocessBx2d(bx, by, layout, i);
            }

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirY]))
            {
                postprocessBy2d(bx, by, layout, i);
            }
        }

        else if constexpr (dimension == 3)
        {
            auto meshSize = layout.meshSize();

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirX]))
            {
                postprocessBx3d(bx, by, bz, meshSize, layout, i);
            }

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirY]))
            {
                postprocessBy3d(bx, by, bz, meshSize, layout, i);
            }

            for (auto const& i : phare_box_from<dimension>(fine_field_box[dirZ]))
            {
                postprocessBz3d(bx, by, bz, meshSize, layout, i);
            }
        }
    }


    static void postprocessBx1d(auto& bx, auto const& layout, core::Point<int, dimension> idx)
    {
        auto locIdx = layout.AMRToLocal(idx);
        auto ix     = locIdx[dirX];
        if (idx[dirX] % 2 != 0)
            bx(ix) = 0.5 * (bx(ix - 1) + bx(ix + 1));
    }

    static void postprocessBx2d(auto& bx, auto& by, auto const& layout,
                                core::Point<int, dimension> idx)
    {
        auto locIdx = layout.AMRToLocal(idx);
        auto ix     = locIdx[dirX];
        auto iy     = locIdx[dirY];
        //                            | <- here with offset = 1
        //                          -- --
        //                            | <- or here with offset = 0
        if (idx[dirX] % 2 != 0)
        {
            // If dual no offset, ie primal for the field we are actually
            // modifying, but dual for the field we are indexing to compute
            // second and third order terms, then the formula reduces to offset
            // = 1
            int xoffset = 1;
            int yoffset = (idx[dirY] % 2 == 0) ? 0 : 1;

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
        auto locIdx = layout.AMRToLocal(idx);
        auto ix     = locIdx[dirX];
        auto iy     = locIdx[dirY];
        //                            |
        //  here with offset = 0 -> -- -- <- or here with offset = 1
        //                            |
        if (idx[dirY] % 2 != 0)
        {
            int xoffset = (idx[dirX] % 2 == 0) ? 0 : 1;
            int yoffset = 1;

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
        auto Dx = meshSize[dirX];
        auto Dy = meshSize[dirY];
        auto Dz = meshSize[dirZ];

        auto locIdx = layout.AMRToLocal(idx);
        auto ix     = locIdx[dirX];
        auto iy     = locIdx[dirY];
        auto iz     = locIdx[dirZ];

        if (idx[dirX] % 2 != 0)
        {
            int xoffset = 1;
            int yoffset = (idx[dirY] % 2 == 0) ? 0 : 1;
            int zoffset = (idx[dirZ] % 2 == 0) ? 0 : 1;

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
        auto Dx = meshSize[dirX];
        auto Dy = meshSize[dirY];
        auto Dz = meshSize[dirZ];

        auto locIdx = layout.AMRToLocal(idx);
        auto ix     = locIdx[dirX];
        auto iy     = locIdx[dirY];
        auto iz     = locIdx[dirZ];

        if (idx[dirY] % 2 != 0)
        {
            int xoffset = (idx[dirX] % 2 == 0) ? 0 : 1;
            int yoffset = 1;
            int zoffset = (idx[dirZ] % 2 == 0) ? 0 : 1;

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
        auto Dx = meshSize[dirX];
        auto Dy = meshSize[dirY];
        auto Dz = meshSize[dirZ];

        auto locIdx = layout.AMRToLocal(idx);
        auto ix     = locIdx[dirX];
        auto iy     = locIdx[dirY];
        auto iz     = locIdx[dirZ];

        if (idx[dirZ] % 2 != 0)
        {
            int xoffset = (idx[dirX] % 2 == 0) ? 0 : 1;
            int yoffset = (idx[dirY] % 2 == 0) ? 0 : 1;
            int zoffset = 1;

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
