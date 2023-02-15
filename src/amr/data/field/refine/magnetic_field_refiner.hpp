#ifndef PHARE_MAGNETIC_FIELD_REFINER_HPP
#define PHARE_MAGNETIC_FIELD_REFINER_HPP

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
class MagneticFieldRefiner
{
public:
    MagneticFieldRefiner(std::array<core::QtyCentering, dimension> const& centering,
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
                if (fineIndex[0] % 2 == 0)
                {
                    fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
                }
                else
                {
                    fineField(locFineIdx[dirX])
                        = 0.5
                          * (coarseField(locCoarseIdx[dirX]) + coarseField(locCoarseIdx[dirX] + 1));
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
                fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
            }
        }




        else if constexpr (dimension == 2)
        {
            if (centerings_[dirX] == core::QtyCentering::primal
                and centerings_[dirY] == core::QtyCentering::dual)
            {
                // Bx
                if (fineIndex[dirX] % 2 == 0)
                {
                    // we're on a coarse X face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
                }
                else
                {
                    // we're on a fine X face, take the average of the coarse value
                    // between the two surrounding faces
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = 0.5
                          * (coarseField(coarseIdx[dirX], coarseIdx[dirY])
                             + coarseField(coarseIdx[dirX] + 1, coarseIdx[dirY]));
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::primal)
            {
                // By
                if (fineIndex[dirY] % 2 == 0)
                {
                    // we're on a coarse Y face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY]);
                }
                else
                {
                    // we're on a fine Y face, take the average of the coarse value
                    // between the two surrounding faces
                    fineField(locFineIdx[dirX], locFineIdx[dirY])
                        = 0.5
                          * (coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY])
                             + coarseField(locCoarseIdx[dirX], locCoarseIdx[dirY] + 1));
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::dual)
            {
                // Bz
                // we're always on a coarse Z face since there is no dual in z
                // all 4 fine Bz take the coarse Z value
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
                if (fineIndex[dirX] % 2 == 0)
                {
                    // we're on a coarse X face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
                else
                {
                    // we're on a fine X face, take the average of the coarse value
                    // between the two surrounding faces
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix + 1, iy, iz));
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::primal
                     and centerings_[dirZ] == core::QtyCentering::dual)
            {
                // By
                if (fineIndex[dirY] % 2 == 0)
                {
                    // we're on a coarse Y face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
                else
                {
                    // we're on a fine Y face, take the average of the coarse value
                    // between the two surrounding faces
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix, iy + 1, iz));
                }
            }
            else if (centerings_[dirX] == core::QtyCentering::dual
                     and centerings_[dirY] == core::QtyCentering::dual
                     and centerings_[dirZ] == core::QtyCentering::primal)
            {
                // Bz
                if (fineIndex[dirZ] % 2 == 0)
                {
                    // we're on a coarse X face
                    // take the coarse face value
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = coarseField(ix, iy, iz);
                }
                else
                {
                    // we're on a fine Z face, take the average of the coarse value
                    // between the two surrounding faces
                    fineField(locFineIdx[dirX], locFineIdx[dirY], locFineIdx[dirZ])
                        = 0.5 * (coarseField(ix, iy, iz) + coarseField(ix, iy, iz + 1));
                }
            }
        }
    }

private:
    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const& centerings_;
};
} // namespace PHARE::amr


#endif // !PHARE_MAGNETIC_FIELD_REFINER_HPP
