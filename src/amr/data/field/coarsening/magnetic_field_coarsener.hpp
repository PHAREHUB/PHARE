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

// clang-format off

// clang-format on

template<std::size_t dim>
template<std::size_t D, typename FieldT>
typename std::enable_if<D == 3, void>::type
MagneticFieldCoarsener<dim>::coarsen(Point_t const fineStartIndex, FieldT const& fineField,
                                     FieldT& coarseField, Point_t const coarseIndex)
{
    auto& ci  = coarseIndex;
    auto& cf  = coarseField;
    auto& fsi = fineStartIndex;
    auto& ff  = fineField;

    using Fn     = std::function<void()>;
    using Fn_arr = std::array<Fn, 2>;
    using A0_arr = std::array<Fn_arr, 2>;

    // X = 0 = primal
    // Y = 1 = dual
    std::array<A0_arr, 2> const fns{
        A0_arr{Fn_arr{[]() { // XXX
                          throw std::runtime_error("no magnetic field should end up here");
                      },
                      [&]() { // XXY
                          double constexpr static V = .5;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) + //
                                    ff(fsi[dirX], fsi[dirY], fsi[dirZ] + 1))
                                   * V;
                      }},
               Fn_arr{[&]() { // XYX
                          double constexpr static V = .5;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) + //
                                    ff(fsi[dirX], fsi[dirY] + 1, fsi[dirZ]))
                                   * V;
                      },
                      [&]() { // XYY
                          double constexpr static V = .5;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) + //
                                    ff(fsi[dirX], fsi[dirY] + 1, fsi[dirZ] + 1))
                                   * V;
                      }}},
        A0_arr{Fn_arr{[&]() {
                          // YXX
                          double constexpr static V = .5;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) + //
                                    ff(fsi[dirX] + 1, fsi[dirY], fsi[dirZ]))
                                   * V;
                      },
                      [&]() {
                          // YXY
                          double constexpr static V = .5;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) + //
                                    ff(fsi[dirX] + 1, fsi[dirY], fsi[dirZ] + 1))
                                   * V;
                      }},
               Fn_arr{[&]() {
                          // YYX
                          double constexpr static V = .5;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) + //
                                    ff(fsi[dirX] + 1, fsi[dirY] + 1, fsi[dirZ]))
                                   * V;
                      },
                      [&]() {
                          // YYY
                          double constexpr static V = .125;

                          cf(ci) = (ff(fsi[dirX], fsi[dirY], fsi[dirZ]) +         //
                                    ff(fsi[dirX] + 1, fsi[dirY], fsi[dirZ]) +     //
                                    ff(fsi[dirX], fsi[dirY] + 1, fsi[dirZ]) +     //
                                    ff(fsi[dirX], fsi[dirY], fsi[dirZ] + 1) +     //
                                    ff(fsi[dirX] + 1, fsi[dirY] + 1, fsi[dirZ]) + //
                                    ff(fsi[dirX] + 1, fsi[dirY], fsi[dirZ] + 1) + //
                                    ff(fsi[dirX], fsi[dirY] + 1, fsi[dirZ] + 1) + //
                                    ff(fsi[dirX] + 1, fsi[dirY] + 1, fsi[dirZ] + 1))
                                   * V;
                      }}}};

    auto cast = [](auto const& qty) {
        return static_cast<std::underlying_type_t<core::QtyCentering>>(qty);
    };
    fns[cast(centering_[0])][cast(centering_[1])][cast(centering_[2])]();
}

} // namespace PHARE::amr
#endif
