#ifndef PHARE_MOMENTS_COARSENER
#define PHARE_MOMENTS_COARSENER

#include "amr/amr_constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"
#include "amr/resources_manager/amr_utils.hpp"


#include <SAMRAI/hier/Box.h>
#include <cstddef>
#include <stdexcept>

namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;
/** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
 * given coarse node
 *
 * A MomentsCoarsener object is created each time the refine() method of the
 * FieldCoarsenOperator is called and its operator() is called for each coarse index.
 * It is the default coarsening policy and used for moments, using primal-primal-primal centering.
 *
 */
template<std::size_t dimension>
class MomentsCoarsener
{
public:
    MomentsCoarsener(std::array<core::QtyCentering, dimension> const centering,
                     SAMRAI::hier::Box const& sourceBox, SAMRAI::hier::Box const& destinationBox,
                     SAMRAI::hier::IntVector const& /*ratio*/)
        : centering_{centering}
        , sourceBox_{sourceBox}
        , destinationBox_{destinationBox}

    {
    }

    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField,
                    core::Point<int, dimension> coarseIndex)
    {
        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        core::Point<int, dimension> fineStartIndex;

        for (auto i = std::size_t{0}; i < dimension; ++i)
        {
            fineStartIndex[i] = coarseIndex[i] * refinementRatio;
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);

        if constexpr (dimension == 1)
        {
            assert(centering_[dirX] == core::QtyCentering::primal);
            coarseField(coarseIndex[dirX]) = fineField(fineStartIndex[dirX]);
        }

        if constexpr (dimension == 2)
        {
            assert(centering_[dirX] == core::QtyCentering::primal
                   and centering_[dirY] == core::QtyCentering::primal);
            coarseField(coarseIndex[dirX], coarseIndex[dirY])
                = fineField(fineStartIndex[dirX], fineStartIndex[dirY]);
        }
        else if constexpr (dimension == 3)
        {
            assert(centering_[dirX] == core::QtyCentering::primal
                   and centering_[dirY] == core::QtyCentering::primal
                   and centering_[dirZ] == core::QtyCentering::primal);
            coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                = fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ]);
        }
    }

private:
    std::array<core::QtyCentering, dimension> const centering_;
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
};

} // namespace PHARE::amr

#endif
