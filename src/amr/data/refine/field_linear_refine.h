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
class UniformIntervalPartitionWeight
{
public:
    UniformIntervalPartitionWeight(QtyCentering centering, int ratio, std::size_t nbrPoints);


    std::vector<double> const& getWeights() const { return weights_; }

private:
    std::vector<double> weights_;
};



template<std::size_t dimension>
class FieldLinearRefineIndexesAndWeights
{
public:
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
        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            if (centering[iDir] == QtyCentering::primal && evenRatio[iDir])
            {
                std::size_t const nbrPoints = static_cast<std::size_t>(ratio(iDir)) + 1;

                UniformIntervalPartitionWeight weight{centering[iDir], ratio(iDir), nbrPoints};

                weights_[iDir] = weight.getWeights();
            }
            else
            {
                std::size_t const nbrPoints = static_cast<std::size_t>(ratio(iDir));

                UniformIntervalPartitionWeight weight{centering[iDir], ratio(iDir), nbrPoints};

                weights_[iDir] = weight.getWeights();
            }
        }

        for (std::size_t iDir = dirX; iDir < dimension; ++iDir)
        {
            if (centering[iDir] == QtyCentering::primal)
            {
                shifts_[iDir] = 0.;
            }
            else
            {
                shifts_[iDir] = 0. - halfRatio[iDir];
            }
        }
    }

    Point<int, dimension> computeStartIndexes(Point<int, dimension> fineIndex) const
    {
        Point<int, dimension> coarseIndex{fineIndex};

        coarseIndex[dirX] = static_cast<int>(static_cast<double>(coarseIndex[dirX] + shifts_[dirX])
                                             / ratio_(dirX));

        if constexpr (dimension > 1)
        {
            coarseIndex[dirY] = static_cast<int>(
                static_cast<double>(coarseIndex[dirY] + shifts_[dirY]) / ratio_(dirY));
        }

        if constexpr (dimension > 2)
        {
            coarseIndex[dirZ] = static_cast<int>(
                static_cast<double>(coarseIndex[dirZ] + shifts_[dirZ]) / ratio_(dirZ));
        }

        return coarseIndex;
    }

    std::array<std::vector<double>, dimension> const& getWeights() const { return weights_; }

    Point<int, dimension> computeIndexesWeight(Point<int, dimension> fineIndex) const
    {
        Point<int, dimension> indexesWeights{fineIndex};

        indexesWeights[dirX] %= ratio_(dirX);

        if constexpr (dimension > 1)
        {
            indexesWeights[dirY] %= ratio_(dirY);
        }
        if constexpr (dimension > 2)
        {
            indexesWeights[dirZ] %= ratio_(dirZ);
        }

        return indexesWeights;
    }

private:
    SAMRAI::hier::IntVector const ratio_;

    std::array<std::vector<double>, dimension> weights_;
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


    template<typename GridLayoutT, typename FieldT>
    void operator()(FieldT const& sourceField, FieldT& destinationField,
                    GridLayoutT const& sourceLayout, GridLayoutT const& destinationLayout,
                    Point<int, dimension> fineIndex)
    {
        TBOX_ASSERT(sourceField.physicalQuantities() == coarseField.physicalQuantities());

        //

        Point<int, dimension> coarseStartIndex = indexesAndWeights_.computeStartIndexes(fineIndex);

        coarseStartIndex = AMRToLocal(coarseStartIndex, coarseBox_);

        Point<int, dimension> iWeight{indexesAndWeights_.computeIndexesWeight(fineIndex)};

        fineIndex = AMRToLocal(fineIndex, fineBox_);

        double fieldWeight = 0.;

        if constexpr (dimension == 1)
        {
            auto const& xStartIndex = coarseStartIndex[dirX];

            auto const& xWeight = weights_[dirX][iWeight[dirX]];

            std::array<double, 2> xWeights{
                {1. - xWeight, xWeight}}; // TODO we want this all the time



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

            auto const& xWeight = weights_[dirX][iWeight[dirX]];
            auto const& yWeight = weights_[dirY][iWeight[dirY]];


            std::array<double, 2> xWeights{
                {1. - xWeight, xWeight}}; // TODO we want this all the time

            std::array<double, 2> yWeights{
                {1. - yWeight, yWeight}}; // TODO we want this all the time

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

            auto const& xWeight = weights_[dirX][iWeight[dirX]];
            auto const& yWeight = weights_[dirY][iWeight[dirY]];
            auto const& zWeight = weights_[dirY][iWeight[dirZ]];


            std::array<double, 2> xWeights{
                {1. - xWeight, xWeight}}; // TODO we want this all the time

            std::array<double, 2> yWeights{
                {1. - yWeight, yWeight}}; // TODO we want this all the time

            std::array<double, 2> zWeights{
                {1. - zWeight, zWeight}}; // TODO we want this all the time

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
    std::array<std::vector<double>, dimension> const& weights_;
};



} // namespace PHARE

#endif
