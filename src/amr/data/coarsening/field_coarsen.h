#ifndef PHARE_FIELD_COARSEN_H
#define PHARE_FIELD_COARSEN_H

#include "data/field/field.h"
#include "data/grid/gridlayoutdefs.h"
#include "hybrid/hybrid_quantities.h"
#include "tools/amr_utils.h"
#include "utilities/constants.h"

#include <SAMRAI/hier/Box.h>

#include <cmath>
#include <vector>



namespace PHARE
{
class Weight
{
public:
    explicit Weight(std::size_t nbrPoints)
    {
        assert(nbrPoints > 1); // we want to have at least 2 points for coarsening operations
        computeWeights_(nbrPoints);
    }

    std::vector<double> getWeights() { return weights_; }

private:
    std::vector<double> weights_;

    double findX_(std::size_t nbrPoints) const;

    void computeWeights_(std::size_t nbrPoints);
};

template<std::size_t dimension>
class IndexesAndWeights
{
public:
    IndexesAndWeights(std::array<QtyCentering, dimension> centering,
                      SAMRAI::hier::IntVector const& ratio)
        : ratio_{ratio}
    {
        std::array<bool, dimension> evenRatio;
        std::array<int, dimension> halfRatio;

        for (int iDir = dirX; iDir < dimension; ++iDir)
        {
            evenRatio[iDir] = ratio(iDir) % 2 == 0;
            halfRatio[iDir] = ratio(iDir) / 2;
        }


        // compute weights for each dimension;
        for (int iDir = dirX; iDir < dimension; ++iDir)
        {
            if (centering[iDir] == QtyCentering::primal || !evenRatio[iDir])
            {
                std::size_t const nbrPoints = 1 + 2 * static_cast<std::size_t>(halfRatio[iDir]);

                Weight weight{nbrPoints};

                weights_[iDir] = weight.getWeights();
            }
            else
            {
                std::size_t const nbrPoints = static_cast<std::size_t>(ratio(iDir));
                Weight weight{nbrPoints};

                weights_[iDir] = weight.getWeights();
            }
        }

        // compute shift for each dimension
        for (int iDir = dirX; iDir < dimension; ++iDir)
        {
            if (centering[iDir] == QtyCentering::primal)
            {
                shifts_[iDir] = 0 - halfRatio[iDir];
            }
            else
            {
                shifts_[iDir] = 0;
            }
        }
    }

    Point<int, dimension> computeStartIndexes(Point<int, dimension> coarseIndex)
    {
        coarseIndex[dirX] = coarseIndex[dirX] * this->ratio_(dirX) + shifts_[dirX];
        if constexpr (dimension > 1)
        {
            coarseIndex[dirY] = coarseIndex[dirY] * this->ratio_(dirY) + shifts_[dirY];
        }
        if constexpr (dimension > 2)
        {
            coarseIndex[dirZ] = coarseIndex[dirZ] * this->ratio_(dirZ) + shifts_[dirZ];
        }

        return coarseIndex;
    }


    std::array<std::vector<double>, dimension> const& getWeights() const { return weights_; }

private:
    SAMRAI::hier::IntVector const ratio_;

    std::array<std::vector<double>, dimension> weights_;
    Point<int, dimension> shifts_;
};


/** @brief Given a dimension, compute the coarsening
 * operation on one coarseIndex, using weights and indexes from
 * the IndexesAndWeights objects
 *
 */
template<std::size_t dimension>
class CoarseField
{
public:
    CoarseField(std::array<QtyCentering, dimension> const& centering,
                SAMRAI::hier::Box const& fineBox, SAMRAI::hier::Box const& coarseBox,
                SAMRAI::hier::IntVector const& ratio)
        : indexesAndWeights_{centering, ratio}
        , fineBox_{fineBox}
        , coarseBox_{coarseBox}
        , weights_{indexesAndWeights_.getWeights()}
    {
    }

    /** @brief apply the coarsening operation of the fineField to the coarseField
     *   at the amr indexes coarseIndex. it is assumed that the fineField will have
     * enought ghost cell for the operation(ghostStencil should be in accordance to
     * the number of ghost for a fineField quantities)
     *
     */
    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField,

                    Point<int, dimension> coarseIndex)
    {
        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantities() == coarseField.physicalQuantities());

        Point<int, dimension> fineStartIndex = indexesAndWeights_.computeStartIndexes(coarseIndex);

        fineStartIndex = AMRToLocal(fineStartIndex, fineBox_);


        coarseIndex = AMRToLocal(coarseIndex, coarseBox_);

        double fieldWeight = 0.;

        if constexpr (dimension == 1)
        {
            auto const& xStartIndex = fineStartIndex[dirX];


            auto const& xWeights = weights_[dirX];


            for (int ix = 0; ix < xWeights.size(); ++ix)
            {
                fieldWeight += fineField(xStartIndex + ix) * xWeights[ix];
            }

            coarseField(coarseIndex[dirX]) = fieldWeight;
        }
        else if constexpr (dimension == 2)
        {
            auto const& xStartIndex = fineStartIndex[dirX];
            auto const& yStartIndex = fineStartIndex[dirY];


            auto const& xWeights = weights_[dirX];
            auto const& yWeights = weights_[dirY];

            for (int ix = 0; ix < xWeights.size(); ++ix)
            {
                double Yinterp = 0.;
                for (int iy = 0; iy < yWeights.size(); ++iy)
                {
                    Yinterp += fineField(xStartIndex + ix, yStartIndex + iy) * yWeights[iy];
                }

                fieldWeight += Yinterp * xWeights[ix];
            }
            coarseField(coarseIndex[dirX], coarseIndex[dirY]) = fieldWeight;
        }
        else if constexpr (dimension == 3)
        {
            auto const& xStartIndex = fineStartIndex[dirX];
            auto const& yStartIndex = fineStartIndex[dirY];
            auto const& zStartIndex = fineStartIndex[dirZ];


            auto const& xWeights = weights_[dirX];
            auto const& yWeights = weights_[dirY];
            auto const& zWeights = weights_[dirZ];



            for (int ix = 0; ix < xWeights.size(); ++ix)
            {
                double Yinterp = 0.;

                for (int iy = 0; iy < yWeights.size(); ++iy)
                {
                    double Zinterp = 0.;
                    for (auto iz = 0; iz < zWeights.size(); ++iz)
                    {
                        Zinterp += fineField(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)
                                   * zWeights[iz];
                    }
                    Yinterp += Zinterp * yWeights[iy];
                }
                fieldWeight += Yinterp * xWeights[ix];
            }

            coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ]) = fieldWeight;
        }
    }



private:
    IndexesAndWeights<dimension> indexesAndWeights_;

    SAMRAI::hier::Box const fineBox_;
    SAMRAI::hier::Box const coarseBox_;

    std::array<std::vector<double>, dimension> const weights_;
};



} // namespace PHARE

#endif
