#ifndef PHARE_ELECTRIC_FIELD_REFINER_HPP
#define PHARE_ELECTRIC_FIELD_REFINER_HPP

#include <SAMRAI/hier/Box.h>

#include "amr/resources_manager/amr_utils.hpp"
#include "core/utilities/constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/constants.hpp"

#include <cstddef>

namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;

template<std::size_t dimension>
class ElectricFieldRefiner
{
public:
    ElectricFieldRefiner(std::array<core::QtyCentering, dimension> const& centering,
                         SAMRAI::hier::Box const& destinationGhostBox,
                         SAMRAI::hier::Box const& sourceGhostBox,
                         SAMRAI::hier::IntVector const& ratio)
        : fineBox_{destinationGhostBox}
        , coarseBox_{sourceGhostBox}
        , centerings_{centering}
    {
    }


    /** @brief Given a sourceField , a destinationField, and a fineIndex compute the
     */
    template<typename FieldT>
    void operator()(FieldT const& coarseField, FieldT& fineField,
                    core::Point<int, dimension> fineIndex)
    {
        TBOX_ASSERT(coarseField.physicalQuantity() == fineField.physicalQuantity());

        auto locFineIdx = AMRToLocal(fineIndex, fineBox_);
        auto coarseIdx{fineIndex};
        for (auto& idx : coarseIdx)
            idx = idx / 2;
        auto locCoarseIdx = AMRToLocal(coarseIdx, coarseBox_);


        if constexpr (dimension == 1)
        {
            // if dual, then Ex
            // if even fine index, we're on top of coarse, we take 100% coarse overlaped fieldValue
            // e.g. fineIndex==100, we take coarse[100/2]
            // e.g. fineIndex==101i  we take coarse[101/2] = coarse[50] as well
            //
            //          49           50             51
            //    o     x     o      x       o      x      o        Ex on coarse
            //    o  x  o  x  o   x  o   x   o   x  o   x  o        Ex on fine
            //       98    99   100     101     102    103
            //
            //
            if (centerings_[0] == core::QtyCentering::dual)
            {
                fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
            }
            // if primal, then Ey,Ez we just copy the coarse value
            // since they are on top of the fine one.
            else
            {
                fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
            }
        }
        else if constexpr (dimension == 2)
        {
            // Ey
            if (centerings_[dirX] == core::QtyCentering::primal
                and centerings_[dirY] == core::QtyCentering::dual)
            {
                if (fineIndex[dirX] % 2 == 0)
                {
                    // we're on a coarse edge at j=1/4 and 3/4
                    // we thus copy the coarse at j=1/2
                    // both fine Ey e.g. at j=100 and 101 will take j=50 on coarse
                    // so no need to look at whether jfine is even or odd
                    // just take the value at the local coarse index
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
                }
                else
                {
                    // we're on a fine edge in between two coarse ones
                    // we take the average
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = 0.5
                          * (coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY])
                             + coarseField(locCoarseIdx[dirX] + 1, locCoarseIdx[dirY]));
                }
            }
            // this is Ex
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::primal)
            {
                if (fineIndex[dirY] % 2 == 0)
                {
                    // we're on a fine edge shared with coarse mesh
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
                }
                else
                {
                    // we're on a fine edge in between two coarse edges
                    // we take the average
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = 0.5
                          * (coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY])
                             + coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY] + 1));
                }
            }
            // and this is now Ez
            else if (centerings_[dirX] == core::QtyCentering::primal
                     and centerings_[dirY] == core::QtyCentering::primal)
            {
                // Ez is always on top of coarse values in 2D
                // so we just copy the coarse value
                fineField(locFineIdx[dirX], locFineIdx[dirY])
                    = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
            }
        }


        else if constexpr (dimension == 3)
        {
            auto ix = locCoarseIdx[dirX];
            auto iy = locCoarseIdx[dirY];
            auto iz = locCoarseIdx[dirZ];

            // Ex
            if (centerings_[dirX] == core::QtyCentering::dual
                and centerings_[dirY] == core::QtyCentering::primal
                and centerings_[dirZ] == core::QtyCentering::primal)
            {
                // here we share an X edge with coarse
                // just copy the coarse value
                if (fineIndex[dirY] % 2 == 0 and fineIndex[dirZ] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
                // we share the Y face but not the Z face
                // we must be one of the 2 X fine edges on a Y face
                // thus we take the average of the two surrounding edges at Z and Z+DZ
                else if (fineIndex[dirY] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix, iy, iz + 1));
                }
                // we share a Z face but not the Y face
                // we must be one of the 2 X fine edges on a Z face
                // we thus take the average of the two X edges at y and y+dy
                else if (fineIndex[dirZ] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix, iy + 1, iz));
                }
                else
                {
                    // we don't share any face thus we're on one of the 2 middle X edges
                    // we take the average of the 4 surrounding X averages
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.25 * (coarseField(ix, iy, iz) + coarseField(ix, iy + 1, iz))
                          + 0.25 * (coarseField(ix, iy, iz + 1) + coarseField(ix, iy + 1, iz + 1));
                }
            }
            // now this is Ey
            else if (centerings_[dirX] == core::QtyCentering::primal
                     and centerings_[dirY] == core::QtyCentering::dual
                     and centerings_[dirZ] == core::QtyCentering::primal)
            {
                // we're on coarse X and coarse Z face, i.e. share a Y coarse edge
                if (fineIndex[dirX] % 2 == 0 and fineIndex[dirZ] % 2 == 0)
                {
                    // we thus just copy the coarse value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
                // now we only have same X face, but not (else) the Z face
                // so we're a new fine Y edge in between two coarse Y edges
                // on a X coarse face. Thus we average the two surrounding Y edges on that X face.
                else if (fineIndex[dirX] % 2 == 0)
                {
                    // we're on a fine Y face, take the average of the coarse value
                    // between the two surrounding faces
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix, iy, iz + 1));
                }
                // we're on a Z coarse face, but not on a X coarse face
                // we thus must be one of the 2 Y edges on a 2 face
                // and thus we take the average of the 2 Y edges at X and X+dX
                else if (fineIndex[dirZ] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix + 1, iy, iz));
                }
                // now we're not on any of the coarse faces
                // so we must be one of the two Y edge in the middle of the cell
                // we thus average over the 4 Y edges of the coarse cell
                else
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.25
                          * (coarseField(ix, iy, iz) + coarseField(ix + 1, iy, iz)
                             + coarseField(ix, iy, iz + 1) + coarseField(ix, iy + 1, iz + 1));
                }
            }
            // now let's do Ez
            else if (centerings_[dirX] == core::QtyCentering::primal
                     and centerings_[dirY] == core::QtyCentering::primal
                     and centerings_[dirZ] == core::QtyCentering::dual)
            {
                // we're on a X and a Y coarse face, we thus share the Z coarse edge
                // we thus copy the coarse value
                if (fineIndex[dirX] % 2 == 0 and fineIndex[dirY] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
                // here we're on a coarse X face, but not a Y face
                // we must be 1 of the 2 Z edges on a X face
                // thus we average the 2 surrounding Z coarse edges at Y and Y+dY
                else if (fineIndex[dirX] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix, iy + 1, iz));
                }
                // here we're on a coarse Y face, but not a X face
                // we must be 1 of the 2 Z edges on a Y face
                // thus we average the 2 surrounding Z coarse edges at X and X+dX
                else if (fineIndex[dirY] % 2 == 0)
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix + 1, iy, iz));
                }
                // we're not on any coarse face thus must be one of the 2 Z edges
                // in the middle of the coarse cell
                // we therefore take the average of the 4 surrounding Z edges
                else
                {
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.25
                          * (coarseField(ix, iy, iz) + coarseField(ix + 1, iy, iz)
                             + coarseField(ix, iy + 1, iz + 1) + coarseField(ix + 1, iy + 1, iz));
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


#endif // !PHARE_ELECTRIC_FIELD_REFINER_HPP
