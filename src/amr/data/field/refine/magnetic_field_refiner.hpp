#ifndef PHARE_MAGNETIC_FIELD_REFINER_HPP
#define PHARE_MAGNETIC_FIELD_REFINER_HPP

#include <SAMRAI/hier/Box.h>

#include "core/utilities/constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"

#include <cstddef>

namespace PHARE::amr
{

template<std::size_t dimension>
class MagneticFieldRefiner
{
public:
    MagneticFieldRefiner(std::array<core::QtyCentering, dimension> const& centering,
                         SAMRAI::hier::Box const& destinationGhostBox,
                         SAMRAI::hier::Box const& sourceGhostBox,
                         SAMRAI::hier::IntVector const& ratio)
        : indexesAndWeights_{centering, ratio}
        , fineBox_{destinationGhostBox}
        , coarseBox_{sourceGhostBox}
    {
    }


    /** @brief Given a sourceField , a destinationField, and a fineIndex compute the
     */
    template<typename FieldT>
    void operator()(FieldT const& sourceField, FieldT& destinationField,
                    core::Point<int, dimension> fineIndex)
    {
        TBOX_ASSERT(sourceField.physicalQuantity() == destinationField.physicalQuantity());

        // First we get the coarseStartIndex for a given fineIndex
        // then we get the index in weights table for a given fineIndex.
        // After that we get the local index of coarseStartIndex and fineIndex.

        // Finally we can compute the interpolation


        core::Point<int, dimension> coarseStartIndex
            = indexesAndWeights_.coarseStartIndex(fineIndex);
        core::Point<int, dimension> iWeight{indexesAndWeights_.computeWeightIndex(fineIndex)};

        coarseStartIndex = AMRToLocal(coarseStartIndex, coarseBox_);
        fineIndex        = AMRToLocal(fineIndex, fineBox_);

        double fieldValue = 0.;




        if constexpr (dimension == 1)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];

            auto const& xWeights         = indexesAndWeights_.weights(core::Direction::X);
            auto const& leftRightWeights = xWeights[iWeight[dirX]];

            for (std::size_t iShiftX = 0; iShiftX < leftRightWeights.size(); ++iShiftX)
            {
                fieldValue += sourceField(xStartIndex + iShiftX) * leftRightWeights[iShiftX];
            }
            destinationField(fineIndex[dirX]) = fieldValue;
        }




        else if constexpr (dimension == 2)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];
            auto const& yStartIndex = coarseStartIndex[dirY];

            auto const& xWeights = indexesAndWeights_.weights(core::Direction::X);
            auto const& yWeights = indexesAndWeights_.weights(core::Direction::Y);

            auto const& xLeftRightWeights = xWeights[iWeight[dirX]];
            auto const& yLeftRightWeights = yWeights[iWeight[dirY]];

            for (std::size_t iShiftX = 0; iShiftX < xLeftRightWeights.size(); ++iShiftX)
            {
                double Yinterp = 0.;
                for (std::size_t iShiftY = 0; iShiftY < yLeftRightWeights.size(); ++iShiftY)
                {
                    Yinterp += sourceField(xStartIndex + iShiftX, yStartIndex + iShiftY)
                               * yLeftRightWeights[iShiftY];
                }
                fieldValue += Yinterp * xLeftRightWeights[iShiftX];
            }

            destinationField(fineIndex[dirX], fineIndex[dirY]) = fieldValue;
        }




        else if constexpr (dimension == 3)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];
            auto const& yStartIndex = coarseStartIndex[dirY];
            auto const& zStartIndex = coarseStartIndex[dirZ];

            auto const& xWeights = indexesAndWeights_.weights(core::Direction::X);
            auto const& yWeights = indexesAndWeights_.weights(core::Direction::Y);
            auto const& zWeights = indexesAndWeights_.weights(core::Direction::Z);

            auto const& xLeftRightWeights = xWeights[iWeight[dirX]];
            auto const& yLeftRightWeights = yWeights[iWeight[dirY]];
            auto const& zLeftRightWeights = zWeights[iWeight[dirZ]];


            for (std::size_t iShiftX = 0; iShiftX < xLeftRightWeights.size(); ++iShiftX)
            {
                double Yinterp = 0.;
                for (std::size_t iShiftY = 0; iShiftY < yLeftRightWeights.size(); ++iShiftY)
                {
                    double Zinterp = 0.;
                    for (std::size_t iShiftZ = 0; iShiftZ < zLeftRightWeights.size(); ++iShiftZ)
                    {
                        Zinterp += sourceField(xStartIndex + iShiftX, yStartIndex + iShiftY,
                                               zStartIndex + iShiftZ)
                                   * zLeftRightWeights[iShiftZ];
                    }
                    Yinterp += Zinterp * yLeftRightWeights[iShiftY];
                }
                fieldValue += Yinterp * xLeftRightWeights[iShiftX];
            }

            destinationField(fineIndex[dirX], fineIndex[dirY], fineIndex[dirZ]) = fieldValue;
        }
    }

private:
    FieldRefineIndexesAndWeights<dimension> const indexesAndWeights_;
    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
};
} // namespace PHARE::amr


#endif // !PHARE_MAGNETIC_FIELD_REFINER_HPP
