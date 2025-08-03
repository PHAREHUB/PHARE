#ifndef PHARE_MHD_FLUX_COARSENER
#define PHARE_MHD_FLUX_COARSENER


#include "core/def/phare_mpi.hpp"

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/constants.hpp"


#include <SAMRAI/hier/Box.h>
#include <cassert>
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
class MHDFluxCoarsener
{
public:
    MHDFluxCoarsener(std::array<core::QtyCentering, dimension> const centering,
                     SAMRAI::hier::Box const& sourceBox, SAMRAI::hier::Box const& destinationBox,
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
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        core::Point<int, dimension> fineStartIndex;

        for (auto i = std::size_t{0}; i < dimension; ++i)
        {
            fineStartIndex[i] = coarseIndex[i] * this->ratio_;
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);

        if constexpr (dimension == 1)
        {
            assert(centering_[dirX] == core::QtyCentering::primal
                   && "MHD flux should be primal in x in 1D");

            coarseField(coarseIndex[dirX]) = fineField(fineStartIndex[dirX]);
        }

        if constexpr (dimension == 2)
        {
            if (centering_[dirX] == core::QtyCentering::primal)
            {
                assert(centering_[dirY] == core::QtyCentering::dual
                       && "MHD flux in x direction should be dual in y");

                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            else if (centering_[dirY] == core::QtyCentering::primal)
            {
                assert(centering_[dirX] == core::QtyCentering::dual
                       && "MHD flux in y direction should be dual in x");

                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            else
            {
                throw std::runtime_error("no MHD flux should end up here");
            }
        }
        else if constexpr (dimension == 3)
        {
            if (centering_[dirX] == core::QtyCentering::primal)
            {
                assert(centering_[dirY] == core::QtyCentering::dual
                       && centering_[dirZ] == core::QtyCentering::dual
                       && "MHD flux in x direction should be dual in y and z");
                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1,
                                     fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY],
                                     fineStartIndex[dirZ] + 1)
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1,
                                     fineStartIndex[dirZ] + 1));
            }
            else if (centering_[dirY] == core::QtyCentering::primal)
            {
                assert(centering_[dirX] == core::QtyCentering::dual
                       && centering_[dirZ] == core::QtyCentering::dual
                       && "MHD flux in y direction should be dual in x and z");
                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY],
                                     fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY],
                                     fineStartIndex[dirZ] + 1)
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY],
                                     fineStartIndex[dirZ] + 1));
            }
            else if (centering_[dirZ] == core::QtyCentering::primal)
            {
                assert(centering_[dirX] == core::QtyCentering::dual
                       && centering_[dirY] == core::QtyCentering::dual
                       && "MHD flux in z direction should be dual in x and y");
                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                    = 0.25
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY],
                                     fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1,
                                     fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1,
                                     fineStartIndex[dirZ]));
            }
            else
            {
                throw std::runtime_error("no MHD flux should end up here");
            }
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
