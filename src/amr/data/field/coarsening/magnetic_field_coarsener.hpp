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
 * WARNING:
 * the following assumes where B is, i.e. Yee layout centering
 * as it only does faces pirmal-dual, dual-primal and dual-dual
 *
 */
template<std::size_t dim>
class MagneticFieldCoarsener
{
    using Point_t = core::Point<int, dim>;

public:
    MagneticFieldCoarsener(std::array<core::QtyCentering, dim> const centering,
                           SAMRAI::hier::Box const& sourceBox,
                           SAMRAI::hier::Box const& destinationBox,
                           SAMRAI::hier::IntVector const& ratio)
        : centering_{centering}
        , sourceBox_{sourceBox}
        , destinationBox_{destinationBox}
    {
    }


    template<typename FieldT>
    void operator()(FieldT const& fineField, FieldT& coarseField, core::Point<int, dim> coarseIndex)
    {
        // For the moment we only take the case of field with the same centering
        TBOX_ASSERT(fineField.physicalQuantity() == coarseField.physicalQuantity());

        coarsen<dim>(fine_start_index(coarseIndex), fineField, coarseField,
                     AMRToLocal(coarseIndex, destinationBox_));
    }

private:
    auto fine_start_index(core::Point<int, dim> const coarseIndex) const
    {
        core::Point<int, dim> fineStartIndex;
        fineStartIndex[dirX] = coarseIndex[dirX] * this->ratio_;
        if constexpr (dim > 1)
        {
            fineStartIndex[dirY] = coarseIndex[dirY] * this->ratio_;
            if constexpr (dim > 2)
            {
                fineStartIndex[dirZ] = coarseIndex[dirZ] * this->ratio_;
            }
        }
        return AMRToLocal(fineStartIndex, sourceBox_);
    }


    template<std::size_t D, typename FieldT>
    typename std::enable_if<D == 1, void>::type
    coarsen(Point_t const fineStartIndex, FieldT const& fineField, FieldT& coarseField,
            Point_t const coarseIndex);

    template<std::size_t D, typename FieldT>
    typename std::enable_if<D == 2, void>::type
    coarsen(Point_t const fineStartIndex, FieldT const& fineField, FieldT& coarseField,
            Point_t const coarseIndex);

    template<std::size_t D, typename FieldT>
    typename std::enable_if<D == 3, void>::type
    coarsen(Point_t const fineStartIndex, FieldT const& fineField, FieldT& coarseField,
            Point_t const coarseIndex);


    std::array<core::QtyCentering, dim> const centering_;
    SAMRAI::hier::Box const sourceBox_;
    SAMRAI::hier::Box const destinationBox_;
    static int constexpr ratio_ = 2;
};

template<std::size_t dim>
template<std::size_t D, typename FieldT>
typename std::enable_if<D == 1, void>::type
MagneticFieldCoarsener<dim>::coarsen(Point_t const fineStartIndex, FieldT const& fineField,
                                     FieldT& coarseField, Point_t const coarseIndex)
{
    // in 1D div(B) is automatically satisfied so using this coarsening
    // opertor is probably not better than the default one, but we do that
    // for some constistency
    // coarse flux is equal to fine flux and we're 1D so tehre is flux partitionned
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

template<std::size_t dim>
template<std::size_t D, typename FieldT>
typename std::enable_if<D == 2, void>::type
MagneticFieldCoarsener<dim>::coarsen(Point_t const fineStartIndex, FieldT const& fineField,
                                     FieldT& coarseField, Point_t const coarseIndex)
{
    //
    //  > == fine bx
    //  >>> = coarse bx
    //
    //   o______________o______________o
    //   |              |              |
    //   |              |              |
    //   >              >              >     0.5
    //   |              |              |     |
    //   |              |              |     |
    //  >>>_____________o_____________>>>  <--
    //   |              |              |     |
    //   |              |              |     |
    //   >              >              >     0.5
    //   |              |              |
    //   |              |              |
    //   o______________o______________o
    //
    // Bx is (primal,dual)
    if (centering_[dirX] == core::QtyCentering::primal
        and centering_[dirY] == core::QtyCentering::dual)
    {
        coarseField(coarseIndex[dirX], coarseIndex[dirY])
            = 0.5
              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1));
    }
    //
    //  > == fine by
    //  >>> = coarse by
    //
    //   o______________o______________o
    //   |              |              |
    //   |              |              |
    //   |              |              |     0.5
    //   |              |              |     |
    //   |              |              |     |
    //   o______________o_____________ o   <--
    //   |              |              |     |
    //   |              |              |     |
    //   |              |              |     0.5
    //   |              |              |
    //   |              ^              |
    //   o_______^______^_______^_______o
    //                  ^
    //          0.5             0.5
    //
    // By is (dual,primal)
    else if (centering_[dirX] == core::QtyCentering::dual
             and centering_[dirY] == core::QtyCentering::primal)
    {
        coarseField(coarseIndex[dirX], coarseIndex[dirY])
            = 0.5
              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY])
                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY]));
    }
    //
    //  fbz == fine bz
    //  CBZ = coarse bz
    //
    //   o______________o______________o
    //   |              |              |
    //   |              |              |
    //   |    fbz       |      fbz     |     0.25
    //   |              |              |     |
    //   |              |              |     |
    //   o_____________CBZ____________ o   <--
    //   |              |              |     |
    //   |              |              |     |
    //   |     fbz      |     fbz      |     0.25
    //   |              |              |
    //   |              |              |
    //   o______________o_______________o
    //
    //          0.25             0.25
    //
    // By is (dual,dual)
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


template<std::size_t dim>
template<std::size_t D, typename FieldT>
typename std::enable_if<D == 3, void>::type
MagneticFieldCoarsener<dim>::coarsen(Point_t const fineStartIndex, FieldT const& fineField,
                                     FieldT& coarseField, Point_t const coarseIndex)
{
    // there are 6 coarse faces on a 3D cell but only 3 correspond to the given coarseIndex

    //  fbx == fine bx
    //  CBX = coarse bx
    //
    //   o______________o______________o
    //   |              |              |
    //   |              |              |
    //   |    fbx       |      fbx     |     0.25
    //   |              |              |     |
    //   |              |              |     |
    //   o_____________CbX____________ o   <--
    //   |              |              |     |
    //   |              |              |     |
    //   |     fbx      |     fbx      |     0.25
    //   |              |              |
    //   |              |              |
    //   o______________o_______________o   ----> z
    //
    //        0.25             0.25
    // Bx is (primal,dual,dual) so average in y and z directions
    if (centering_[dirX] == core::QtyCentering::primal
        and centering_[dirY] == core::QtyCentering::dual
        and centering_[dirZ] == core::QtyCentering::dual)
    {
        coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
            = 0.25
              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1, fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ] + 1)
                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1,
                             fineStartIndex[dirZ] + 1));
    }

    //  fby == fine by
    //  CBy = coarse by
    //
    //   o______________o______________o
    //   |              |              |
    //   |              |              |
    //   |    fby       |      fby     |     0.25
    //   |              |              |     |
    //   |              |              |     |
    //   o_____________CBY____________ o   <--
    //   |              |              |     |
    //   |              |              |     |
    //   |     fby      |     fby      |     0.25
    //   |              |              |
    //   |              |              |
    //   o______________o_______________o   ----> x
    //
    //        0.25             0.25
    // By is (dual,primal,dual) so average in x and z directions
    else if (centering_[dirX] == core::QtyCentering::dual
             and centering_[dirY] == core::QtyCentering::primal
             and centering_[dirZ] == core::QtyCentering::dual)
    {
        coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
            = 0.25
              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY], fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ] + 1)
                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY],
                             fineStartIndex[dirZ] + 1));
    }

    //  fbz == fine bz
    //  CBz = coarse bz
    //
    //   o______________o______________o
    //   |              |              |
    //   |              |              |
    //   |    fbz       |      fbz     |     0.25
    //   |              |              |     |
    //   |              |              |     |
    //   o_____________CBZ____________ o   <--
    //   |              |              |     |
    //   |              |              |     |
    //   |     fbz      |     fbz      |     0.25
    //   |              |              |
    //   |              |              |
    //   o______________o______________o   ----> x
    //
    //        0.25             0.25
    // Bz is (dual,dual,primal) so average in x and y directions
    else if (centering_[dirX] == core::QtyCentering::dual
             and centering_[dirY] == core::QtyCentering::dual
             and centering_[dirZ] == core::QtyCentering::primal)
    {
        coarseField(coarseIndex[dirX], coarseIndex[dirY], coarseIndex[dirZ])
            = 0.25
              * (fineField(fineStartIndex[dirX], fineStartIndex[dirY], fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY], fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX], fineStartIndex[dirY] + 1, fineStartIndex[dirZ])
                 + fineField(fineStartIndex[dirX] + 1, fineStartIndex[dirY] + 1,
                             fineStartIndex[dirZ]));
    }
    else
    {
        throw std::runtime_error("no magnetic field should end up here");
    }
}

} // namespace PHARE::amr
#endif
