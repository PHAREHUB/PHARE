#ifndef PHARE_MAGNETIC_FIELD_REGRIDER_HPP
#define PHARE_MAGNETIC_FIELD_REGRIDER_HPP


#include "core/def/phare_mpi.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include <SAMRAI/hier/Box.h>

#include <cmath>
#include <cstddef>

namespace PHARE::amr
{

/** \brief Refines the magnetic components from a coarse mesh to fine faces shared with the coarse
 * ones.
 *
 * This refinement operator works for magnetic field components dispatched following the Yee layout.
 * It sets the values of fine components only on faces shared with coarse faces.
 * The fine faces values are set equal to that of the coarse shared one (order 0 interpolation).
 * inner fine faces are set by the MagneticRefinePatchStrategy
 */
template<std::size_t dimension>
class MagneticFieldRegrider
{
public:
    MagneticFieldRegrider(std::array<core::QtyCentering, dimension> const& centering,
                          SAMRAI::hier::Box const& destinationGhostBox,
                          SAMRAI::hier::Box const& sourceGhostBox,
                          SAMRAI::hier::IntVector const& /*ratio*/)
        : fineBox_{destinationGhostBox}
        , coarseBox_{sourceGhostBox}
        , centerings_{centering}
    {
    }


    // magnetic field refinement is made so to conserve the divergence of B
    // it simply copies the value of the magnetic field existing on a coarse face
    // onto the 2 (1D), 4 (2/3D) colocated fine faces. This way the total flux on
    // these fine faces equals that on the overlaped coarse face.
    // see fujimoto et al. 2011 :  doi:10.1016/j.jcp.2011.08.002
    template<typename FieldT>
    void operator()(FieldT const& coarseField, FieldT& fineField,
                    core::Point<int, dimension> fineIndex)
    {
        TBOX_ASSERT(coarseField.physicalQuantity() == fineField.physicalQuantity());

        using core::dirX;
        using core::dirY;
        using core::dirZ;

        auto const locFineIdx   = AMRToLocal(fineIndex, fineBox_);
        auto const coarseIdx    = toCoarseIndex(fineIndex);
        auto const locCoarseIdx = AMRToLocal(coarseIdx, coarseBox_);


        if constexpr (dimension == 1)
        {
            // if primal, i.e. Bx :
            // if even fine index, we're on top of coarse, we take 100% coarse overlaped fieldValue
            // e.g. fineIndex==100, we take coarse[100/2]
            // if odd fine index, we take 50% of surrounding coarse nodes
            // e.g. fineIndex == 101, we take 0.5(coarse(101/2)+coarse(101/2+1))
            //
            //   49          50              51            52
            //    o           o              o             o      Bx on coarse
            //    x     x     x      x       o      x      x      Bx on fine
            //   98     99   100    101     102    103    104
            //
            //
            if (centerings_[0] == core::QtyCentering::primal)
            {
                if (fineIndex[0] % 2 == 0 && std::isnan(fineField(locFineIdx[dirX])))
                {
                    fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
                }
            }
            // dual case, By, Bz
            //          49           50             51
            //    o     +     o      +       o      +       o      Byz on coarse : +
            //    o  +  o  +  o   +  o   +   o   +  o   +   o      Byz on fine   : +
            //      98    99     100    101     102    103
            //
            // 100 takes 50 = 100/2
            // 101 takes 50 = 101/2
            else
            {
                if (std::isnan(fineField(locFineIdx[dirX])))
                    fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
            }
        }




        else if constexpr (dimension == 2)
        {
            if (centerings_[dirX] == core::QtyCentering::primal
                and centerings_[dirY] == core::QtyCentering::dual)
            {
                // Bx
                if (fineIndex[dirX] % 2 == 0
                    && std::isnan(fineField(locFineIdx[dirX], locFineIdx[dirY])))
                {
                    // we're on a coarse X face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::primal)
            {
                // By
                if (fineIndex[dirY] % 2 == 0
                    && std::isnan(fineField(locFineIdx[dirX], locFineIdx[dirY])))
                {
                    // we're on a coarse Y face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::dual)
            {
                // Bz
                // we're always on a coarse Z face since there is no dual in z
                // all 4 fine Bz take the coarse Z value
                if (std::isnan(fineField(locFineIdx[dirX], locFineIdx[dirY])))
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
            }
        }


        else if constexpr (dimension == 3)
        {
            auto ix = locCoarseIdx[dirX];
            auto iy = locCoarseIdx[dirY];
            auto iz = locCoarseIdx[dirZ];

            if (centerings_[dirX] == core::QtyCentering::primal
                and centerings_[dirY] == core::QtyCentering::dual
                and centerings_[dirZ] == core::QtyCentering::dual)
            {
                // Bx
                if (fineIndex[dirX] % 2 == 0
                    && std::isnan(fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])))
                {
                    // we're on a coarse X face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::primal
                     and centerings_[dirZ] == core::QtyCentering::dual)
            {
                // By
                if (fineIndex[dirY] % 2 == 0
                    && std::isnan(fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])))
                {
                    // we're on a coarse Y face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::dual
                     and centerings_[dirZ] == core::QtyCentering::primal)
            {
                // Bz
                if (fineIndex[dirZ] % 2 == 0
                    && std::isnan(fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])))
                {
                    // we're on a coarse X face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
            }
        }
    }

private:
    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
};
} // namespace PHARE::amr


#endif // !PHARE_MAGNETIC_FIELD_REFINER_HPP
