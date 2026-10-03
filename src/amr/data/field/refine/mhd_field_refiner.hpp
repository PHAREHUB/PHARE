#ifndef PHARE_MHD_FIELD_REFINER_HPP
#define PHARE_MHD_FIELD_REFINER_HPP


#include "phare_mpi.hpp"

#include <SAMRAI/hier/Box.h>

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/refine/any_field_refiner.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"

#include <cstddef>

namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;

template<std::size_t dimension>
class MHDFieldRefiner : public AnyFieldRefiner<MHDFieldRefiner<dimension>, dimension>
{
public:
    using Super = AnyFieldRefiner<MHDFieldRefiner<dimension>, dimension>;
    using Super::Super;

protected:
    using Super::centerings_;
    using Super::coarseBox_;
    using Super::fineBox_;

public:
    // electric field refinement strategy follows
    // fujimoto et al. 2011 :  doi:10.1016/j.jcp.2011.08.002
    template<typename FieldT>
    void refine(FieldT const& coarseField, FieldT& fineField, auto const fineIndex)
    {
        TBOX_ASSERT(coarseField.physicalQuantity() == fineField.physicalQuantity());

        auto const locFineIdx   = AMRToLocal(fineIndex, fineBox_);
        auto const coarseIdx    = toCoarseIndex(fineIndex);
        auto const locCoarseIdx = AMRToLocal(coarseIdx, coarseBox_);

        refine_(coarseField, fineField, locFineIdx, locCoarseIdx);
    }

private:
    template<typename FieldT>
        requires(dimension == 1)
    void refine_(FieldT const& coarseField, FieldT& fineField, auto const& locFineIdx,
                 auto const& locCoarseIdx)
    {
        assert(centerings_[dirX] == core::QtyCentering::dual
               && "MHD field should be primal in x in 1D");

        if (auto& fineVal = fineField(locFineIdx); std::isnan(fineVal))
            fineVal = coarseField(locCoarseIdx);
    }

    template<typename FieldT>
        requires(dimension == 2)
    void refine_(FieldT const& coarseField, FieldT& fineField, auto const& locFineIdx,
                 auto const& locCoarseIdx)
    {
        assert(centerings_[dirX] == core::QtyCentering::dual
               && centerings_[dirY] == core::QtyCentering::dual
               && "MHD field should be dual in x and y in 2D");

        if (auto& fineVal = fineField(locFineIdx); std::isnan(fineVal))
            fineVal = coarseField(locCoarseIdx);
    }


    template<typename FieldT>
        requires(dimension == 3)
    void refine_(FieldT const& coarseField, FieldT& fineField, auto const& locFineIdx,
                 auto const& locCoarseIdx)
    {
        assert(centerings_[dirX] == core::QtyCentering::dual
               && centerings_[dirY] == core::QtyCentering::dual
               && centerings_[dirZ] == core::QtyCentering::dual
               && "MHD field should be dual in x, y and z in 3D");

        if (auto& fineVal = fineField(locFineIdx); std::isnan(fineVal))
            fineVal = coarseField(locCoarseIdx);
    }
};
} // namespace PHARE::amr


#endif // PHARE_MHD_FIELD_REFINER_HPP
