#ifndef PHARE_FIELD_MOMENTS_REFINER_HPP
#define PHARE_FIELD_MOMENTS_REFINER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/field/field.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include "field_linear_refine.hpp"

#include <SAMRAI/hier/Box.h>

#include <array>
#include <vector>


namespace PHARE
{
namespace amr
{
    /**@brief a FieldRefiner is an object that is used to get the value of a field at a fine AMR
     * index from coarse data
     *
     * The FieldRefiner is created each time a refinement is needed by the FieldRefinementOperator
     * and its operator() is used for each fine index onto which we want to get the value from the
     * coarse field.
     */
    template<std::size_t dimension>
    class FieldMomentsRefiner
    {
    public:
        FieldMomentsRefiner(std::array<core::QtyCentering, dimension> const& centering,
                            SAMRAI::hier::Box const& destinationGhostBox,
                            SAMRAI::hier::Box const& sourceGhostBox,
                            SAMRAI::hier::IntVector const& ratio)
            : indexesAndWeights_{centering, ratio}
            , fineBox_{destinationGhostBox}
            , coarseBox_{sourceGhostBox}
        {
        }


        /** @brief Given a sourceField , a destinationField, and a fineIndex compute the
         * interpolation from the coarseField(sourceField) to the fineFiled(destinationField) at the
         * fineIndex index
         *
         *
         * Strategy :
         * - for a given fineIndex, we first compute the associated CoarseIndex
         * - the two coarse indexes to get coarse values are then coarseIndex and coarseIndex+1
         * - the weights are pre-computed by the FieldRefineIndexesAndWeights object
         * - we just have to know which one to use, depending on where the fineIndex is in the
         * coarse cell
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
} // namespace amr
} // namespace PHARE


#endif
