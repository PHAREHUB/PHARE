#ifndef PHARE_ELECTRIC_FIELD_COARSENER
#define PHARE_ELECTRIC_FIELD_COARSENER

#include "amr/amr_constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/constants.hpp"
#include "amr/resources_manager/amr_utils.hpp"


#include <SAMRAI/hier/Box.h>
#include <cstddef>
#include <stdexcept>

namespace PHARE::amr
{
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
class ElectricFieldCoarsener
{
public:
    ElectricFieldCoarsener(std::array<core::QtyCentering, dimension> const centering,
                           SAMRAI::hier::Box const& sourceBox,
                           SAMRAI::hier::Box const& destinationBox,
                           SAMRAI::hier::IntVector const& /*ratio*/)
        : centering_{centering}
        , sourceBox_{sourceBox}
        , destinationBox_{destinationBox}

    {
    }

    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField,
                    core::Point<int, dimension> coarseIndex)
    {
        using core::dirX;
        using core::dirY;
        using core::dirZ;

        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        core::Point<int, dimension> fineStartIndex;

        for (auto i = std::size_t{0}; i < dimension; ++i)
        {
            fineStartIndex[i] = coarseIndex[i] * refinementRatio;
        }

        fineStartIndex = AMRToLocal(fineStartIndex, sourceBox_);
        coarseIndex    = AMRToLocal(coarseIndex, destinationBox_);

        if constexpr (dimension == 1)
        {
            if (centering_[dirX] == core::QtyCentering::dual) // ex
            {
                coarseField(coarseIndex[dirX])
                    = 0.5 * (fineField(fineStartIndex[dirX] + 1) + fineField(fineStartIndex[dirX]));
            }
            else if (centering_[dirX] == core::QtyCentering::primal) // ey, ez
            {
                coarseField(coarseIndex[dirX]) = fineField(fineStartIndex[dirX]);
            }
        }

        if constexpr (dimension == 2)
        {
            if (centering_[dirX] == core::QtyCentering::dual
                and centering_[dirY] == core::QtyCentering::primal) // ex
            {
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
            }
            else if (centering_[dirX] == core::QtyCentering::primal
                     and centering_[dirY] == core::QtyCentering::dual) // ey
            {
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
            }
            else if (centering_[dirX] == core::QtyCentering::primal
                     and centering_[dirY] == core::QtyCentering::primal) // ez
            {
                coarseField(coarseIndex[dirX], coarseIndex[dirY])
                    = fineField(fineStartIndex[dirX], fineStartIndex[dirY]);
            }
            else
            {
                throw std::runtime_error("no electric field should end up here");
            }
        }
        else if constexpr (dimension == 3)
        {
            if (centering_[dirX] == core::QtyCentering::dual
                and centering_[dirY] == core::QtyCentering::primal
                and centering_[dirZ] == core::QtyCentering::primal) // ex
            {
                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY],
                                     fineStartIndex[dirZ]));
            }
            else if (centering_[dirX] == core::QtyCentering::primal
                     and centering_[dirY] == core::QtyCentering::dual
                     and centering_[dirZ] == core::QtyCentering::primal) // ey
            {
                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1,
                                     fineStartIndex[dirZ]));
            }
            else if (centering_[dirX] == core::QtyCentering::primal
                     and centering_[dirY] == core::QtyCentering::primal
                     and centering_[dirZ] == core::QtyCentering::dual) // ez
            {
                coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
                    = 0.5
                      * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                         + fineField(fineStartIndex[dirX], fineStartIndex[dirY],
                                     fineStartIndex[dirZ] + 1));
            }
            else
            {
                throw std::runtime_error("no electric field should end up here");
            }
        }
    }

private:
    std::array<core::QtyCentering, dimension> const centering_;
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
};

} // namespace PHARE::amr

#endif
