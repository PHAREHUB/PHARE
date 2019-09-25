#ifndef PHARE_FIELD_COARSEN_H
#define PHARE_FIELD_COARSEN_H

#include "coarsen_weighter.h"
#include "data/field/field.h"
#include "data/grid/gridlayoutdefs.h"
#include "hybrid/hybrid_quantities.h"
#include "resources_manager/amr_utils.h"
#include "utilities/constants.h"



#include <SAMRAI/hier/Box.h>

#include <cmath>
#include <vector>



namespace PHARE
{
namespace amr
{
    template<std::size_t dimension>
    class FieldCoarsenIndexesAndWeights
    {
    private:
        std::array<std::size_t, dimension>
        nbrFinePoints_(std::array<core::QtyCentering, dimension> centering,
                       SAMRAI::hier::IntVector const& ratio) const
        {
            std::array<std::size_t, dimension> nbrPoints;

            for (auto iDir = 0u; iDir < dimension; ++iDir)
            {
                if (centering[iDir] == core::QtyCentering::primal && ratio(iDir) % 2 == 0)
                {
                    nbrPoints[iDir] = static_cast<std::size_t>(ratio(iDir)) + 1;
                }
                else
                {
                    nbrPoints[iDir] = static_cast<std::size_t>(ratio(iDir));
                }
            }
            return nbrPoints;
        }


    public:
        FieldCoarsenIndexesAndWeights(std::array<core::QtyCentering, dimension> centering,
                                      SAMRAI::hier::IntVector const& ratio)
            : ratio_{ratio}
            , weighters_{make_weighters(nbrFinePoints_(centering, ratio),
                                        std::make_index_sequence<dimension>{})}
        {
            std::array<bool, dimension> evenRatio;
            std::array<int, dimension> halfRatio;

            for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
            {
                evenRatio[iDir] = ratio(iDir) % 2 == 0;
                halfRatio[iDir] = ratio(iDir) / 2;
            }

            // In primal centering, coarseIndex*RF is RF/2 on the right of
            // the startFineIndex, so we store -Rf/2 in shift to later add to cearseIndex*RF
            // to get the startIndex
            for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
            {
                if (centering[iDir] == core::QtyCentering::primal)
                {
                    shifts_[iDir] = 0 - halfRatio[iDir];
                }
                else
                {
                    shifts_[iDir] = 0;
                }
            }
        }

        core::Point<int, dimension>
        computeStartIndexes(core::Point<int, dimension> const& coarseIndex)
        {
            core::Point<int, dimension> fineIndex{coarseIndex};
            fineIndex[dirX] = coarseIndex[dirX] * this->ratio_(dirX) + shifts_[dirX];
            if constexpr (dimension > 1)
            {
                fineIndex[dirY] = coarseIndex[dirY] * this->ratio_(dirY) + shifts_[dirY];
            }
            if constexpr (dimension > 2)
            {
                fineIndex[dirZ] = coarseIndex[dirZ] * this->ratio_(dirZ) + shifts_[dirZ];
            }

            return fineIndex;
        }


        std::vector<double> const& weights(core::Direction dir) const
        {
            return weighters_[static_cast<std::size_t>(dir)].weights();
        }

    private:
        SAMRAI::hier::IntVector const ratio_;
        std::array<CoarsenWeighter, dimension> weighters_;
        core::Point<int, dimension> shifts_;
    };

} // namespace amr


} // namespace PHARE

#endif
