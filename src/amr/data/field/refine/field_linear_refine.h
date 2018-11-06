#ifndef PHARE_FIELD_LINEAR_REFINE_H
#define PHARE_FIELD_LINEAR_REFINE_H


#include "data/field/field.h"
#include "data/grid/gridlayoutdefs.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"

#include <SAMRAI/hier/Box.h>

#include <array>
#include <unordered_map>
#include <vector>


namespace PHARE
{
/** @brief  This class will contain uniform spaced distance in the interval 0,1
 * the distance is from the left index.
 */
class UniformIntervalPartitionWeight
{
public:
    UniformIntervalPartitionWeight(QtyCentering centering, std::size_t ratio,
                                   std::size_t nbrPoints);


    std::vector<double> const& getUniformDistances() const { return distances_; }

private:
    std::vector<double> distances_;
    std::vector<int> relativeDualIndexes;
};




template<std::size_t dimension>
class FieldRefineIndexesAndWeights
{
public:
    /** @brief Given a centering in each directions and a ratio, initialize weights and shifts
     *  for later use. (the vector of weights will be the same regardless of the fineIndex)
     * it is which index of the weights that will be used depends on the fineIndex, and
     * also which coarseIndex to start for refine operation
     *
     */
    FieldRefineIndexesAndWeights(std::array<QtyCentering, dimension> centering,
                                 SAMRAI::hier::IntVector const& ratio)
        : ratio_{ratio}
    {
        std::array<bool, dimension> isEvenRatio;

        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            isEvenRatio[iDir] = ratio(iDir) % 2 == 0;
        }

        // compute weights for each directions
        // number of points is ratio + 1 if we are primal or odd ratio
        // and ratio otherwise
        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            // here we extract the distances of the left index
            // and then compute the weights : 1.-distance for the left one
            // and distance for the right index
            // The number of points depends on the centering, for primal or odd ratio
            // it is ratio + 1 , for dual with evenRatio it is ratio
            auto nbrPoints = static_cast<std::size_t>(ratio(iDir));

            UniformIntervalPartitionWeight distances{
                centering[iDir], static_cast<std::size_t>(ratio(iDir)), nbrPoints};

            weights_[iDir].reserve(distances.getUniformDistances().size());
            for (auto const& distance : distances.getUniformDistances())
            {
                weights_[iDir].emplace_back(std::array<double, 2>{{1. - distance, distance}});
            }
        }

        // this shift will be use to determine which coarseIndexe we take
        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            if (centering[iDir] == QtyCentering::primal)
            {
                shifts_[iDir] = 0.;
            }
            else
            {
                // in case we are dual, we need to shift our fine index of - halfRatio
                // so that after truncating to integer (the index/ratio), we get the correct
                // coarseStartIndex
                shifts_[iDir] = 0.5;
            }
        }
    }




    Point<int, dimension> computeStartIndexes(Point<int, dimension> fineIndex) const
    {
        Point<int, dimension> coarseIndex{fineIndex};

        // here we perform the floating point division, and then we truncate to integer
        coarseIndex[dirX] = static_cast<int>(
            static_cast<double>(fineIndex[dirX] + shifts_[dirX]) / ratio_(dirX) - shifts_[dirX]);

        if constexpr (dimension > 1)
        {
            coarseIndex[dirY] = static_cast<int>(
                static_cast<double>(fineIndex[dirY] + shifts_[dirY]) / ratio_(dirY)
                - shifts_[dirY]);
        }

        if constexpr (dimension > 2)
        {
            coarseIndex[dirZ] = static_cast<int>(
                static_cast<double>(fineIndex[dirZ] + shifts_[dirZ]) / ratio_(dirZ)
                - shifts_[dirZ]);
        }

        return coarseIndex;
    }




    std::array<std::vector<std::array<double, 2>>, dimension> const& getWeights() const
    {
        return weights_;
    }




    /** @brief Compute the index of weigths for a given fineIndex
     *
     */
    Point<int, dimension> computeWeightIndex(Point<int, dimension> fineIndex) const
    {
        Point<int, dimension> indexesWeights{fineIndex};

        indexesWeights[dirX] %= ratio_[dirX];

        if constexpr (dimension > 1)
        {
            indexesWeights[dirY] %= ratio_[dirY];
        }
        if constexpr (dimension > 2)
        {
            indexesWeights[dirZ] %= ratio_[dirZ];
        }

        return indexesWeights;
    }

private:
    SAMRAI::hier::IntVector const ratio_;

    std::array<std::vector<std::array<double, 2>>, dimension> weights_;
    Point<double, dimension> shifts_;
};




template<std::size_t dimension>
class FieldLinearRefine
{
public:
    FieldLinearRefine(std::array<QtyCentering, dimension> const& centering,
                      SAMRAI::hier::Box const& destinationGhostBox,
                      SAMRAI::hier::Box const& sourceGhostBox, SAMRAI::hier::IntVector const& ratio)
        : indexesAndWeights_{centering, ratio}
        , fineBox_{destinationGhostBox}
        , coarseBox_{sourceGhostBox}
        , weights_{indexesAndWeights_.getWeights()}
    {
    }


    /** @brief Given a sourceField , a destinationField, and a fineIndex compute the interpolation
     * from the coarseField(sourceField) to the fineFiled(destinationField) at the fineIndex index
     */
    template<typename FieldT>
    void operator()(FieldT const& sourceField, FieldT& destinationField,
                    Point<int, dimension> fineIndex)
    {
        TBOX_ASSERT(sourceField.physicalQuantities() == coarseField.physicalQuantities());

        // First we get the coarseStartIndex for a given fineIndex
        // then we get the index in weights table for a given fineIndex.
        // After that we get the local index of coarseStartIndex and fineIndex.

        // Finally we can compute the interpolation


        Point<int, dimension> coarseStartIndex = indexesAndWeights_.computeStartIndexes(fineIndex);
        Point<int, dimension> iWeight{indexesAndWeights_.computeWeightIndex(fineIndex)};



        coarseStartIndex = AMRToLocal(coarseStartIndex, coarseBox_);
        fineIndex        = AMRToLocal(fineIndex, fineBox_);

        double fieldValue = 0.;

        if constexpr (dimension == 1)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];

            auto const& xWeights = weights_[dirX][iWeight[dirX]];


            for (std::size_t iShiftX = 0; iShiftX < xWeights.size(); ++iShiftX)
            {
                fieldValue += sourceField(xStartIndex + iShiftX) * xWeights[iShiftX];
            }


            destinationField(fineIndex[dirX])
                = fieldValue; // TODO : field should take a MeshIndex/Point (choose and kill the
                              // other?)
        }
        else if constexpr (dimension == 2)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];
            auto const& yStartIndex = coarseStartIndex[dirY];

            auto const& xWeights = weights_[dirX][iWeight[dirX]];
            auto const& yWeights = weights_[dirY][iWeight[dirY]];



            for (std::size_t iShiftX = 0; iShiftX < xWeights.size(); ++iShiftX)
            {
                double Yinterp = 0.;
                for (std::size_t iShiftY = 0; iShiftY < yWeights.size(); ++iShiftY)
                {
                    Yinterp += sourceField(xStartIndex + iShiftX, yStartIndex + iShiftY)
                               * yWeights[iShiftY];
                }
                fieldValue += Yinterp * xWeights[iShiftX];
            }

            destinationField(fineIndex[dirX], fineIndex[dirY]) = fieldValue;
        }
        else if constexpr (dimension == 3)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];
            auto const& yStartIndex = coarseStartIndex[dirY];
            auto const& zStartIndex = coarseStartIndex[dirZ];

            auto const& xWeights = weights_[dirX][iWeight[dirX]];
            auto const& yWeights = weights_[dirY][iWeight[dirY]];
            auto const& zWeights = weights_[dirY][iWeight[dirZ]];


            for (std::size_t iShiftX = 0; iShiftX < xWeights.size(); ++iShiftX)
            {
                double Yinterp = 0.;
                for (std::size_t iShiftY = 0; iShiftY < yWeights.size(); ++iShiftY)
                {
                    double Zinterp = 0.;
                    for (std::size_t iShiftZ = 0; iShiftZ < zWeights.size(); ++iShiftZ)
                    {
                        Zinterp += sourceField(xStartIndex + iShiftX, yStartIndex + iShiftY,
                                               zStartIndex + iShiftZ)
                                   * zWeights[iShiftZ];
                    }
                    Yinterp += Zinterp * yWeights[iShiftY];
                }
                fieldValue += Yinterp * xWeights[iShiftX];
            }

            destinationField(fineIndex[dirX], fineIndex[dirY], fineIndex[dirZ]) = fieldValue;
        }
    }

private:
    FieldRefineIndexesAndWeights<dimension> const indexesAndWeights_;
    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<std::vector<std::array<double, 2>>, dimension> const& weights_;
};



} // namespace PHARE

#endif
