#ifndef PHARE_MAGNETIC_FIELD_COARSENER
#define PHARE_MAGNETIC_FIELD_COARSENER


#include "core/def/phare_mpi.hpp"

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/constants.hpp"


#include <SAMRAI/hier/Box.h>
#include <stdexcept>

namespace PHARE::amr
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
 *
 *
 * This coarsening operation is defined so to conserve the magnetic flux.
 * This is done by assigning to a magnetic field component on a coarse face, the average
 * of the enclosed fine faces
 *
 */
template<std::size_t dimension>
class MagneticFieldCoarsener
{
public:
    MagneticFieldCoarsener(std::array<core::QtyCentering, dimension> const centering,
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

        fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;

        if constexpr (dimension > 1)
        {
            fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
            if constexpr (dimension > 2)
            {
                fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
            }
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);

        // the following kinda assumes where B is, i.e. Yee layout centering
        // as it only does faces pirmal-dual, dual-primal and dual-dual

        if constexpr (dimension == 1)
        {
            // in 1D div(B) is automatically satisfied so using this coarsening
            // opertor is probably not better than the default one, but we do that
            // for a kind of consistency...
            // coarse flux is equal to fine flux and we're 1D so there is flux partitioned
            // only for By and Bz, Bx is equal to the fine value

            if (centering_[dirX] == core::QtyCentering::primal) // bx
            {
                coarseField(coarseIndex[dirX]) = fineField(fineStartIndex[dirX]);
            }
            else if (centering_[dirX] == core::QtyCentering::dual) // by and bz
            {
                coarseField(coarseIndex[dirX])
                    = 0.5 * (fineField(fineStartIndex[dirX] + 1) + fineField(fineStartIndex[dirX]));
            }
        }

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
            else
            {
                throw std::runtime_error("no magnetic field should end up here");
            }
        }
        else if constexpr (dimension == 3)
        {
            throw std::runtime_error("Not Implemented yet");
        }
    }

private:
    std::array<core::QtyCentering, dimension> const centering_;
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
    static int constexpr ratio_ = 2;
};
} // namespace PHARE::amr
#endif
