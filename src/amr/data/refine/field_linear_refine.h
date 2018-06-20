#ifndef PHARE_FIELD_LINEAR_REFINE_H
#define PHARE_FIELD_LINEAR_REFINE_H


#include "data/field/field.h"
#include "data/grid/gridlayoutdefs.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"

#include <SAMRAI/hier/Box.h>

#include <array>
#include <vector>


namespace PHARE
{
/** @brief  This class will contain uniform spaced distance in the interval 0,1
 * the distance is from the left index.
 *
 */
class UniformIntervalPartitionWeight
{
public:
    UniformIntervalPartitionWeight(QtyCentering centering, int ratio, std::size_t nbrPoints);


    std::vector<double> const& getUniformDistances() const { return distances_; }

private:
    std::vector<double> distances_;
};


template<std::size_t dimension>
class FieldLinearRefineIndexesAndWeights
{
public:
    /** @brief Given a centering in each directions and a ratio, initialize weights and shifts
     *  for later use. (the vector of weights will be the same regardless of the fineIndex)
     * it is which index of the weights that will be used depends on the fineIndex, and
     * also which coarseIndex to start for refine operation
     *
     */
    FieldLinearRefineIndexesAndWeights(std::array<QtyCentering, dimension> centering,
                                       SAMRAI::hier::IntVector const& ratio)
        : ratio_{ratio}
    {
        std::array<bool, dimension> evenRatio;
        std::array<double, dimension> halfRatio;

        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            evenRatio[iDir] = ratio(iDir) % 2 == 0;
            halfRatio[iDir] = ratio(iDir) / 2.;
        }

        // compute weights for each directions
        // number of points is ratio + 1 if we are primal or odd ratio
        // and                 ratio otherwise
        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            // here we extract the distances of the left index
            // and then compute the weights : 1.-distance for the left one
            // and distance for the right index
            // The number of points depends on the centering, for primal or odd ratio
            // it is ratio + 1 , for dual with evenRatio it is ratio
            if (centering[iDir] == QtyCentering::primal || !evenRatio[iDir])
            {
                std::size_t const nbrPoints = static_cast<std::size_t>(ratio(iDir)) + 1;

                UniformIntervalPartitionWeight distances{centering[iDir], ratio(iDir), nbrPoints};

                weights_[iDir].reserve(distances.getUniformDistances().size());
                for (auto const& distance : distances.getUniformDistances())
                {
                    weights_[iDir].emplace_back(std::array<double, 2>{{1. - distance, distance}});
                }
            }
            else
            {
                std::size_t const nbrPoints = static_cast<std::size_t>(ratio(iDir));

                UniformIntervalPartitionWeight distances{centering[iDir], ratio(iDir), nbrPoints};

                weights_[iDir].reserve(distances.getUniformDistances().size());
                for (auto const& distance : distances.getUniformDistances())
                {
                    weights_[iDir].emplace_back(std::array<double, 2>{{1. - distance, distance}});
                }
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
                shifts_[iDir] = 0. - halfRatio[iDir];
            }
        }
    }

    Point<int, dimension> computeStartIndexes(Point<int, dimension> fineIndex) const
    {
        Point<int, dimension> coarseIndex{fineIndex};

        // here we perform the floating point division, and then we truncate to integer
        coarseIndex[dirX]
            = static_cast<int>(static_cast<double>(fineIndex[dirX] + shifts_[dirX]) / ratio_(dirX));

        if constexpr (dimension > 1)
        {
            coarseIndex[dirY] = static_cast<int>(
                static_cast<double>(fineIndex[dirY] + shifts_[dirY]) / ratio_(dirY));
        }

        if constexpr (dimension > 2)
        {
            coarseIndex[dirZ] = static_cast<int>(
                static_cast<double>(fineIndex[dirZ] + shifts_[dirZ]) / ratio_(dirZ));
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

        // Finnaly we can compute the interpolation


        Point<int, dimension> coarseStartIndex = indexesAndWeights_.computeStartIndexes(fineIndex);
        Point<int, dimension> iWeight{indexesAndWeights_.computeWeightIndex(fineIndex)};



        coarseStartIndex = AMRToLocal(coarseStartIndex, coarseBox_);
        fineIndex        = AMRToLocal(fineIndex, fineBox_);

        double fieldWeight = 0.;

        if constexpr (dimension == 1)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];

            auto const& xWeights = weights_[dirX][iWeight[dirX]];


            for (std::size_t ix = 0; ix < xWeights.size(); ++ix)
            {
                fieldWeight += sourceField(xStartIndex + ix) * xWeights[ix];
            }


            destinationField(fineIndex[dirX]) = fieldWeight;
        }
        else if constexpr (dimension == 2)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];
            auto const& yStartIndex = coarseStartIndex[dirY];

            auto const& xWeights = weights_[dirX][iWeight[dirX]];
            auto const& yWeights = weights_[dirY][iWeight[dirY]];



            for (std::size_t ix = 0; ix < xWeights.size(); ++ix)
            {
                double Yinterp = 0.;
                for (std::size_t iy = 0; iy < yWeights.size(); ++iy)
                {
                    Yinterp += sourceField(xStartIndex + ix, yStartIndex + iy) * yWeights[iy];
                }
                fieldWeight += Yinterp * xWeights[ix];
            }

            destinationField(fineIndex[dirX], fineIndex[dirY]) = fieldWeight;
        }
        else if constexpr (dimension == 3)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];
            auto const& yStartIndex = coarseStartIndex[dirY];
            auto const& zStartIndex = coarseStartIndex[dirZ];

            auto const& xWeights = weights_[dirX][iWeight[dirX]];
            auto const& yWeights = weights_[dirY][iWeight[dirY]];
            auto const& zWeights = weights_[dirY][iWeight[dirZ]];


            for (std::size_t ix = 0; ix < xWeights.size(); ++ix)
            {
                double Yinterp = 0.;
                for (std::size_t iy = 0; iy < yWeights.size(); ++iy)
                {
                    double Zinterp = 0.;
                    for (std::size_t iz = 0; iz < zWeights.size(); ++iz)
                    {
                        Zinterp += sourceField(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)
                                   * zWeights[iz];
                    }
                    Yinterp += Zinterp * yWeights[iy];
                }
                fieldWeight += Yinterp * xWeights[ix];
            }

            destinationField(fineIndex[dirX], fineIndex[dirY], fineIndex[dirZ]) = fieldWeight;
        }
    }

private:
    FieldLinearRefineIndexesAndWeights<dimension> const indexesAndWeights_;
    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;
    std::array<std::vector<std::array<double, 2>>, dimension> const& weights_;
};



} // namespace PHARE

#endif
