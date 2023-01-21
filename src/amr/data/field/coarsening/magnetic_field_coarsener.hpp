#ifndef PHARE_MAGNETIC_FIELD_COARSENER
#define PHARE_MAGNETIC_FIELD_COARSENER

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/constants.hpp"

#include <SAMRAI/hier/Box.h>

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
     * A MagneticFieldCoarsener object is created each time the refine() method of the
     * FieldCoarsenOperator is called and its operator() is called for each coarse index.
     * It is the default coarsening policy and used for any field that does not come with
     * specific constraints (such as conserving some property in the coarsening process).
     */
    template<std::size_t dimension>
    class MagneticFieldCoarsener
    {
    public:
        MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const& centering,
                               SAMRAI::hier::Box const& sourceBox,
                               SAMRAI::hier::Box const& destinationBox,
                               SAMRAI::hier::IntVector const& ratio)
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
            if constexpr (dimension == 1)
            {
                fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
            }
            else if constexpr (dimension > 1)
            {
                fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
                if constexpr (dimension > 2)
                {
                    fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
                }
            }

            fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
            coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);


            double coarseValue = 0.;

            // we're not supposed to know B has this specific centering here
            // hard coded for now but should instead used the layout to ask for centering

            if constexpr (dimension == 2)
            {
                if (centering_[dirX] == core::QtyCentering::primal
                    and centering_[dirY] == core::QtyCentering::dual)
                {
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
                }
                else if (centering_[dirX] == core::QtyCentering::dual
                         and centering_[dirY] == core::QtyCentering::primal)
                {
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.5
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
                }
                else if (centering_[dirX] == core::QtyCentering::dual
                         and centering_[dirY] == core::QtyCentering::dual)
                {
                    coarseField(coarseIndex[dirX], coarseIndex[dirY])
                        = 0.25
                          * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY])
                             + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1)
                             + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1));
                }
            }
        }
        std::array<core::QtyCentering, dimension> const& centering_;
        SAMRAI::hier::Box const sourceBox_;
        SAMRAI::hier::Box const destinationBox_;
        static int constexpr ratio_ = 2;
    };
} // namespace amr
} // namespace PHARE
#endif
