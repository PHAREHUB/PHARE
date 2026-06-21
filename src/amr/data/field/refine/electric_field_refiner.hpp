#ifndef PHARE_ELECTRIC_FIELD_REFINER_HPP
#define PHARE_ELECTRIC_FIELD_REFINER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include <SAMRAI/hier/Box.h>

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


    // electric field refinement strategy follows
    // fujimoto et al. 2011 :  doi:10.1016/j.jcp.2011.08.002
    void operator()(auto const& coarseField, auto& fineField, auto const& fineIndex,
                    auto const& coarseIndex, auto const& locFineIdx, auto const& locCoarseIdx,
                    auto& fineVal, auto const& coarseVal)
    {
        TBOX_ASSERT(coarseField.physicalQuantity() == fineField.physicalQuantity());

        if (not std::isnan(fineVal))
            return; // KEEP!

        if constexpr (dimension == 1)
            fineVal = coarseVal;
        else if constexpr (dimension == 2)
            refine2D_(coarseField, fineField, fineIndex, coarseIndex, locFineIdx, locCoarseIdx,
                      fineVal, coarseVal);
        else if constexpr (dimension == 3)
            refine3D_(coarseField, fineField, fineIndex, coarseIndex, locFineIdx, locCoarseIdx,
                      fineVal, coarseVal);
    }

private:
    // knowing we have a refinement ratio of 2, every fine face that has an even index
    // is on top of a coarse face, and every fine face that has an odd index is in between
    // two coarse faces.
    bool onCoarseXFace_(core::Point<int, dimension> const& fineIndex)
    {
        return fineIndex[dirX] % 2 == 0;
    }
    bool onCoarseYFace_(core::Point<int, dimension> const& fineIndex)
    {
        return fineIndex[dirY] % 2 == 0;
    }
    bool onCoarseZFace_(core::Point<int, dimension> const& fineIndex)
    {
        return fineIndex[dirZ] % 2 == 0;
    }


    void refine2D_(auto const& coarseField, auto& fineField, auto const& fineIndex,
                   auto const& coarseIndex, auto const& locFineIdx, auto const& locCoarseIdx,
                   auto& fineVal, auto const& coarseVal)
    {
        // ilc: index local coarse
        auto const ilcx = locCoarseIdx[dirX];
        auto const ilcy = locCoarseIdx[dirY];

        // this is Ex
        if (centerings_[dirX] == core::QtyCentering::dual
            and centerings_[dirY] == core::QtyCentering::primal)
        {
            if (onCoarseYFace_(fineIndex))
            {
                // we're on a fine edge shared with coarse mesh
                // take the coarse face value
                fineVal = coarseVal;
            }
            else
            {
                // we're on a fine edge in between two coarse edges
                // we take the average

                fineVal = 0.5 * (coarseVal + coarseField(ilcx, ilcy + 1));
            }
        }
        // Ey
        else if (centerings_[dirX] == core::QtyCentering::primal
                 and centerings_[dirY] == core::QtyCentering::dual)
        {
            if (onCoarseXFace_(fineIndex))
            {
                // we're on a coarse edge at j=1/4 and 3/4
                // we thus copy the coarse at j=1/2
                // both fine Ey e.g. at j=100 and 101 will take j=50 on coarse
                // so no need to look at whether jfine is even or odd
                // just take the value at the local coarse index

                fineVal = coarseVal;
            }
            else
            {
                // we're on a fine edge in between two coarse ones
                // we take the average

                fineVal = 0.5 * (coarseVal + coarseField(ilcx + 1, ilcy));
            }
        }
        // and this is now Ez
        else if (centerings_[dirX] == core::QtyCentering::primal
                 and centerings_[dirY] == core::QtyCentering::primal)
        {
            if (onCoarseXFace_(fineIndex) and onCoarseYFace_(fineIndex))
            {
                fineVal = coarseVal;
            }
            else if (onCoarseXFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx, ilcy + 1));
            }
            else if (onCoarseYFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx + 1, ilcy));
            }
            else
            {
                fineVal = 0.25
                          * (coarseVal + coarseField(ilcx + 1, ilcy) + coarseField(ilcx, ilcy + 1)
                             + coarseField(ilcx + 1, ilcy + 1));
            }
        }
    }


    void refine3D_(auto const& coarseField, auto& fineField, auto const& fineIndex,
                   auto const& coarseIndex, auto const& locFineIdx, auto const& locCoarseIdx,
                   auto& fineVal, auto const& coarseVal)
    {
        // ilc: index local coarse
        auto const ilcx = locCoarseIdx[dirX];
        auto const ilcy = locCoarseIdx[dirY];
        auto const ilcz = locCoarseIdx[dirZ];

        // Ex
        if (centerings_[dirX] == core::QtyCentering::dual
            and centerings_[dirY] == core::QtyCentering::primal
            and centerings_[dirZ] == core::QtyCentering::primal)
        {
            // here we share an X edge with coarse
            // just copy the coarse value
            if (onCoarseYFace_(fineIndex) and onCoarseZFace_(fineIndex))
            {
                fineVal = coarseVal;
            }
            // we share the Y face but not the Z face
            // we must be one of the 2 X fine edges on a Y face
            // thus we take the average of the two surrounding edges at Z and Z+DZ
            else if (onCoarseYFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx, ilcy, ilcz + 1));
            }
            // we share a Z face but not the Y face
            // we must be one of the 2 X fine edges on a Z face
            // we thus take the average of the two X edges at y and y+dy
            else if (onCoarseZFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx, ilcy + 1, ilcz));
            }
            else
            {
                // we don't share any face thus we're on one of the 2 middle X edges
                // we take the average of the 4 surrounding X averages

                fineVal = 0.25 * (coarseVal + coarseField(ilcx, ilcy + 1, ilcz))
                          + 0.25
                                * (coarseField(ilcx, ilcy, ilcz + 1)
                                   + coarseField(ilcx, ilcy + 1, ilcz + 1));
            }
        }
        // now this is Ey
        else if (centerings_[dirX] == core::QtyCentering::primal
                 and centerings_[dirY] == core::QtyCentering::dual
                 and centerings_[dirZ] == core::QtyCentering::primal)
        {
            // we're on coarse X and coarse Z face, i.e. share a Y coarse edge
            if (onCoarseXFace_(fineIndex) and onCoarseZFace_(fineIndex))
            {
                // we thus just copy the coarse value

                fineVal = coarseVal;
            }
            // now we only have same X face, but not (else) the Z face
            // so we're a new fine Y edge in between two coarse Y edges
            // on a X coarse face. Thus we average the two surrounding Y edges on that X
            // face.
            else if (onCoarseXFace_(fineIndex))
            {
                // we're on a fine X face, but not on a Z face
                // this means we are on a Y edge that lies in between 2 coarse edges
                // at z and z+dz
                // take the average of these 2 coarse value

                fineVal = 0.5 * (coarseVal + coarseField(ilcx, ilcy, ilcz + 1));
            }
            // we're on a Z coarse face, but not on a X coarse face
            // we thus must be one of the 2 Y edges on a Z face
            // and thus we take the average of the 2 Y edges at X and X+dX
            else if (onCoarseZFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx + 1, ilcy, ilcz));
            }
            // now we're not on any of the coarse faces
            // so we must be one of the two Y edge in the middle of the cell
            // we thus average over the 4 Y edges of the coarse cell
            else
            {
                fineVal = 0.25
                          * (coarseVal + coarseField(ilcx + 1, ilcy, ilcz)
                             + coarseField(ilcx, ilcy, ilcz + 1)
                             + coarseField(ilcx + 1, ilcy, ilcz + 1));
            }
        }
        // now let's do Ez
        else if (centerings_[dirX] == core::QtyCentering::primal
                 and centerings_[dirY] == core::QtyCentering::primal
                 and centerings_[dirZ] == core::QtyCentering::dual)
        {
            // we're on a X and a Y coarse face, we thus share the Z coarse edge
            // we thus copy the coarse value
            if (onCoarseXFace_(fineIndex) and onCoarseYFace_(fineIndex))
            {
                fineVal = coarseVal;
            }
            // here we're on a coarse X face, but not a Y face
            // we must be 1 of the 2 Z edges on a X face
            // thus we average the 2 surrounding Z coarse edges at Y and Y+dY
            else if (onCoarseXFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx, ilcy + 1, ilcz));
            }
            // here we're on a coarse Y face, but not a X face
            // we must be 1 of the 2 Z edges on a Y face
            // thus we average the 2 surrounding Z coarse edges at X and X+dX
            else if (onCoarseYFace_(fineIndex))
            {
                fineVal = 0.5 * (coarseVal + coarseField(ilcx + 1, ilcy, ilcz));
            }
            // we're not on any coarse face thus must be one of the 2 Z edges
            // in the middle of the coarse cell
            // we therefore take the average of the 4 surrounding Z edges
            else
            {
                fineVal = 0.25
                          * (coarseVal + coarseField(ilcx + 1, ilcy, ilcz)
                             + coarseField(ilcx, ilcy + 1, ilcz)
                             + coarseField(ilcx + 1, ilcy + 1, ilcz));
            }
        }
    }

    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
};
} // namespace PHARE::amr


#endif // PHARE_ELECTRIC_FIELD_REFINER_HPP
