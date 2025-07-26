#ifndef PHARE_MHD_FLUX_REFINER_HPP
#define PHARE_MHD_FLUX_REFINER_HPP


#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/Box.h>

#include "amr/resources_manager/amr_utils.hpp"
#include "core/utilities/constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"

#include <cstddef>

namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;

template<std::size_t dimension>
class MHDFluxRefiner
{
public:
    MHDFluxRefiner(std::array<core::QtyCentering, dimension> const& centering,
                   SAMRAI::hier::Box const& destinationGhostBox,
                   SAMRAI::hier::Box const& sourceGhostBox, SAMRAI::hier::IntVector const& ratio)
        : fineBox_{destinationGhostBox}
        , coarseBox_{sourceGhostBox}
        , centerings_{centering}
    {
    }


    // electric field refinement strategy follows
    // fujimoto et al. 2011 :  doi:10.1016/j.jcp.2011.08.002
    template<typename FieldT>
    void operator()(FieldT const& coarseField, FieldT& fineField,
                    core::Point<int, dimension> fineIndex)
    {
        TBOX_ASSERT(coarseField.physicalQuantity() == fineField.physicalQuantity());

        auto const locFineIdx   = AMRToLocal(fineIndex, fineBox_);
        auto const coarseIdx    = toCoarseIndex(fineIndex);
        auto const locCoarseIdx = AMRToLocal(coarseIdx, coarseBox_);

        if constexpr (dimension == 1)
            refine1D_(coarseField, fineField, locFineIdx, locCoarseIdx);
        else if constexpr (dimension == 2)
            refine2D_(coarseField, fineField, fineIndex, locFineIdx, locCoarseIdx);
        else if constexpr (dimension == 3)
            refine3D_(coarseField, fineField, fineIndex, locFineIdx, locCoarseIdx);
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


    template<typename FieldT>
    void refine1D_(FieldT const& coarseField, FieldT& fineField,
                   core::Point<int, dimension> const& locFineIdx,
                   core::Point<int, dimension> const& locCoarseIdx)
    {
        assert(centerings_[dirX] == core::QtyCentering::primal
               && "MHD flux should be primal in x in 1D");
        fineField(locFineIdx[dirX]) = coarseField(locCoarseIdx[dirX]);
    }

    template<typename FieldT>
    void refine2D_(FieldT const& coarseField, FieldT& fineField,
                   core::Point<int, dimension> const& fineIndex,
                   core::Point<int, dimension> const& locFineIdx,
                   core::Point<int, dimension> const& locCoarseIdx)
    {
        // ilc: index local coarse
        // ilf: index local fine
        auto const ilcx = locCoarseIdx[dirX];
        auto const ilcy = locCoarseIdx[dirY];
        auto const ilfx = locFineIdx[dirX];
        auto const ilfy = locFineIdx[dirY];


        if (centerings_[dirX] == core::QtyCentering::primal)
        {
            assert(centerings_[dirY] == core::QtyCentering::dual
                   && "MHD flux in x direction should be dual in y");
            if (onCoarseXFace_(fineIndex))
            {
                fineField(ilfx, ilfy) = coarseField(ilcx, ilcy);
            }
            else
            {
                fineField(ilfx, ilfy)
                    = 0.5 * (coarseField(ilcx, ilcy) + coarseField(ilcx + 1, ilcy));
            }
        }
        else if (centerings_[dirY] == core::QtyCentering::primal)
        {
            assert(centerings_[dirX] == core::QtyCentering::dual
                   && "MHD flux in y direction should be dual in x");
            if (onCoarseYFace_(fineIndex))
            {
                fineField(ilfx, ilfy) = coarseField(ilcx, ilcy);
            }
            else
            {
                fineField(ilfx, ilfy)
                    = 0.5 * (coarseField(ilcx, ilcy) + coarseField(ilcx, ilcy + 1));
            }
        }
        else
        {
            throw std::runtime_error(
                "MHDFluxRefiner: no MHD flux should only be x or y centered in 2D");
        }
    }


    template<typename FieldT>
    void refine3D_(FieldT const& coarseField, FieldT& fineField,
                   core::Point<int, dimension> const& fineIndex,
                   core::Point<int, dimension> const& locFineIdx,
                   core::Point<int, dimension> const& locCoarseIdx)
    {
        // ilc: index local coarse
        // ilf: index local fine
        auto const ilcx = locCoarseIdx[dirX];
        auto const ilcy = locCoarseIdx[dirY];
        auto const ilcz = locCoarseIdx[dirZ];
        auto const ilfx = locFineIdx[dirX];
        auto const ilfy = locFineIdx[dirY];
        auto const ilfz = locFineIdx[dirZ];

        if (centerings_[dirX] == core::QtyCentering::primal)
        {
            assert(centerings_[dirY] == core::QtyCentering::dual
                   && centerings_[dirZ] == core::QtyCentering::dual
                   && "MHD flux in x direction should be dual in y and z");
            if (onCoarseXFace_(fineIndex))
            {
                fineField(ilfx, ilfy, ilfz) = coarseField(ilcx, ilcy, ilcz);
            }
            else
            {
                fineField(ilfx, ilfy, ilfz)
                    = 0.5 * (coarseField(ilcx, ilcy, ilcz) + coarseField(ilcx + 1, ilcy, ilcz));
            }
        }
        else if (centerings_[dirY] == core::QtyCentering::primal)
        {
            assert(centerings_[dirX] == core::QtyCentering::dual
                   && centerings_[dirZ] == core::QtyCentering::dual
                   && "MHD flux in y direction should be dual in x and z");
            if (onCoarseYFace_(fineIndex))
            {
                fineField(ilfx, ilfy, ilfz) = coarseField(ilcx, ilcy, ilcz);
            }
            else
            {
                fineField(ilfx, ilfy, ilfz)
                    = 0.5 * (coarseField(ilcx, ilcy, ilcz) + coarseField(ilcx, ilcy + 1, ilcz));
            }
        }
        else if (centerings_[dirZ] == core::QtyCentering::primal)
        {
            assert(centerings_[dirX] == core::QtyCentering::dual
                   && centerings_[dirY] == core::QtyCentering::dual
                   && "MHD flux in z direction should be dual in x and y");
            if (onCoarseZFace_(fineIndex))
            {
                fineField(ilfx, ilfy, ilfz) = coarseField(ilcx, ilcy, ilcz);
            }
            else
            {
                fineField(ilfx, ilfy, ilfz)
                    = 0.5 * (coarseField(ilcx, ilcy, ilcz) + coarseField(ilcx, ilcy, ilcz + 1));
            }
        }
    }

    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
};
} // namespace PHARE::amr


#endif // PHARE_ELECTRIC_FIELD_REFINER_HPP
