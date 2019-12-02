#ifndef PHARE_FIELD_COARSENER_H
#define PHARE_FIELD_COARSENER_H

#include "core/data/grid/gridlayoutdefs.h"
#include "amr/data/field/coarsening/field_coarsen_index_weight.h"
#include "core/utilities/constants.h"

#include <SAMRAI/hier/Box.h>

#include <cstddef>



namespace PHARE
{
namespace amr
{
    using core::dirX;
    using core::dirY;
    using core::dirZ;
    /** @brief This class gives an operator() that performs the coarsening of N fine nodes onto a
     * given coarse node
     *
     * A FieldCoarsener object is created each time the refine() method of the FieldCoarsenOperator
     * is called and its operator() is called for each coarse index.
     */
    template<std::size_t dimension>
    class FieldCoarsener
    {
    public:
        FieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                       SAMRAI::hier::Box const& sourceBox, SAMRAI::hier::Box const& destinationBox,
                       SAMRAI::hier::IntVector const& ratio)
            : indexesAndWeights_{centering, ratio}
            , sourceBox_{sourceBox}
            , destinationBox_{destinationBox}
        {
        }

        /** @brief apply the coarsening operation of the fineField to the coarseField
         *   at the amr indexes coarseIndex. it is assumed that the fineField will have
         * enought ghost cell for the operation(ghostStencil should be in accordance to
         * the number of ghost for a fineField quantities.
         */
        template<typename FieldT>
        void operator()(FieldT const& fineField, FieldT& coarseField,

                        core::Point<int, dimension> coarseIndex)
        {
            // For the moment we only take the case of field with the same centering
            TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

            core::Point<int, dimension> fineStartIndex
                = indexesAndWeights_.computeStartIndexes(coarseIndex);

            fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
            coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


            double coarseValue = 0.;




            if constexpr (dimension == 1)
            {
                auto const& xStartIndex = fineStartIndex[dirX];
                auto const& xWeights    = indexesAndWeights_.weights(core::Direction::X);


                for (std::size_t iShiftX = 0; iShiftX < xWeights.size(); ++iShiftX)
                {
                    coarseValue += fineField(xStartIndex + iShiftX) * xWeights[iShiftX];
                }

                coarseField(coarseIndex[dirX]) = coarseValue;
            }




            else if constexpr (dimension == 2)
            {
                auto const& xStartIndex = fineStartIndex[dirX];
                auto const& yStartIndex = fineStartIndex[dirY];


                auto const& xWeights = indexesAndWeights_.weights(core::Direction::X);
                auto const& yWeights = indexesAndWeights_.weights(core::Direction::Y);

                for (std::size_t iShiftX = 0; iShiftX < xWeights.size(); ++iShiftX)
                {
                    double Yinterp = 0.;
                    for (std::size_t iShiftY = 0; iShiftY < yWeights.size(); ++iShiftY)
                    {
                        Yinterp += fineField(xStartIndex + iShiftX, yStartIndex + iShiftY)
                                   * yWeights[iShiftY];
                    }

                    coarseValue += Yinterp * xWeights[iShiftX];
                }
                coarseField(coarseIndex[dirX], coarseIndex[dirY]) = coarseValue;
            }



            else if constexpr (dimension == 3)
            {
                auto const& xStartIndex = fineStartIndex[dirX];
                auto const& yStartIndex = fineStartIndex[dirY];
                auto const& zStartIndex = fineStartIndex[dirZ];

                auto const& xWeights = indexesAndWeights_.weights(core::Direction::X);
                auto const& yWeights = indexesAndWeights_.weights(core::Direction::Y);
                auto const& zWeights = indexesAndWeights_.weights(core::Direction::Z);

                for (std::size_t iShiftX = 0; iShiftX < xWeights.size(); ++iShiftX)
                {
                    double Yinterp = 0.;

                    for (std::size_t iShiftY = 0; iShiftY < yWeights.size(); ++iShiftY)
                    {
                        double Zinterp = 0.;
                        for (std::size_t iShiftZ = 0; iShiftZ < zWeights.size(); ++iShiftZ)
                        {
                            Zinterp += fineField(xStartIndex + iShiftX, yStartIndex + iShiftY,
                                                 zStartIndex + iShiftZ)
                                       * zWeights[iShiftZ];
                        }
                        Yinterp += Zinterp * yWeights[iShiftY];
                    }
                    coarseValue += Yinterp * xWeights[iShiftX];
                }

                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ]) = coarseValue;
            }
        }



    private:
        //! precompute the indexes and weights to use to coarsen fine values onto a coarse node
        FieldCoarsenIndexesAndWeights<dimension> indexesAndWeights_;
        SAMRAI::hier::Box const sourceBox_;
        SAMRAI::hier::Box const destinationBox_;
    };
} // namespace amr
} // namespace PHARE


#endif
