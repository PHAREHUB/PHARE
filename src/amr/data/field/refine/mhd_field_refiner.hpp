#ifndef PHARE_MHD_FIELD_REFINER_HPP
#define PHARE_MHD_FIELD_REFINER_HPP


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
class MHDFieldRefiner
{
public:
    MHDFieldRefiner(std::array<core::QtyCentering, dimension> const& centering,
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
    template<typename FieldT>
    void refine1D_(FieldT const& coarseField, FieldT& fineField,
                   core::Point<int, dimension> const& locFineIdx,
                   core::Point<int, dimension> const& locCoarseIdx)
    {
        assert(centerings_[dirX] == core::QtyCentering::dual
               && "MHD field should be primal in x in 1D");

        if (std::isnan(fineField(locFineIdx[dirX])))
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

        assert(centerings_[dirX] == core::QtyCentering::dual
               && centerings_[dirY] == core::QtyCentering::dual
               && "MHD field should be dual in x and y in 2D");

        if (std::isnan(fineField(ilfx, ilfy)))
            fineField(ilfx, ilfy) = coarseField(ilcx, ilcy);
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

        assert(centerings_[dirX] == core::QtyCentering::dual
               && centerings_[dirY] == core::QtyCentering::dual
               && centerings_[dirZ] == core::QtyCentering::dual
               && "MHD field should be dual in x, y and z in 3D");

        if (std::isnan(fineField(ilfx, ilfy, ilfz)))
            fineField(ilfx, ilfy, ilfz) = coarseField(ilcx, ilcy, ilcz);
    }

    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<core::QtyCentering, dimension> const centerings_;
};
} // namespace PHARE::amr


#endif // PHARE_ELECTRIC_FIELD_REFINER_HPP
