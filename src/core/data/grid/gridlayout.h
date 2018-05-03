#ifndef PHARE_CORE_DATA_GRID_GRIDLAYOUT_H
#define PHARE_CORE_DATA_GRID_GRIDLAYOUT_H

#include <array>
#include <memory>
#include <vector>


#include "data/field/field.h"
#include "gridlayoutimpl.h"
#include "gridlayoutimplyee.h"

#include "hybrid/hybrid_quantities.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"
#include "utilities/types.h"


namespace PHARE
{
constexpr int centering2int(QtyCentering c)
{
    return static_cast<int>(c);
}

/**
 * @brief Gridlayout is a class used to handle all operations related to
 * a specific grid layout, for instance the Yee layout.
 *
 * A GridLayout is a central object in PHARE. All grid operation need to manipulate
 * a GridLayout. It is used to get start and end indices of ghost nodes, or of
 * nodes of the physical domain. GridLayout is the object that knows how to
 * calculate derivatives, and therefore enable other classes to be written
 * in a very general way, without knowing anything that is specific to the layout
 * of quantities on the grid.
 *
 * The GridLayout class follows the Pimpl idiom to be seen, from client code, as
 * a simple object living in its scope. The implementation pointer is abstract
 * and concrete child classes implement concret layouts, such as Yee.
 */
template<Layout layout, std::size_t dim>
class GridLayout
{
    // ------------------------------------------------------------------------
    //                              PRIVATE
    // ------------------------------------------------------------------------

private:
    uint32 nbDims_{dim}; // dimensionality (1, 2, 3)

    std::array<double, dim> dl_;

    std::array<double, dim> odxyz_;

    std::array<uint32, dim> nbrCell_;

    Point<double, dim> origin_; // origin of the grid
    uint32 interpOrder_;

    std::string layoutName_;

    GridLayoutImpl<layout, dim> impl_; // private implementation

    using error = std::runtime_error; // TODO: find a better error handling
    static const std::string errorInverseMesh;


    // test the validity of the GridLayout construction
    void throwNotValid() const;


    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------

public:
    static const uint32 minNbrCells = 1; // minimum nbr of cells in a
                                         // non-invariant direction

    GridLayout(std::array<double, dim> const& dx, std::array<uint32, dim> const& nbrCells,
               std::string const& layoutName, Point<double, dim> const& orgin,
               uint32 ghostParameter);

    GridLayout(GridLayout const& source) = delete;
    GridLayout(GridLayout&& source)      = delete;

    GridLayout& operator=(GridLayout const& source) = delete;
    GridLayout& operator=(GridLayout&& source) = delete;


    void init(std::array<double, dim> const& dx, std::array<uint32, dim> const& nbrCells,
              Point<double, dim> const& origin);

    Point<double, dim> origin() const { return origin_; }


    std::array<double, dim> dxdydz() const { return dl_; }

    std::array<double, dim> odxyz() const { return odxyz_; }




    std::array<uint32, dim> nbrCellxyz() const { return nbrCell_; }


    uint32 nbDimensions() const { return nbDims_; }

    uint32 order() const { return interpOrder_; }

    std::string layoutName() const { return layoutName_; }

    constexpr static std::array<QtyCentering, dim>
    centering(HybridQuantity::Scalar const& hybridQuantity)
    {
        return GridLayoutImpl<layout, dim>::centering(hybridQuantity);
    }

    uint32 physicalStartIndex(QtyCentering centering, Direction direction) const;
    uint32 physicalStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                              Direction direction) const;
    template<typename NdArrayImpl>
    uint32 physicalStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                              Direction direction) const;

    uint32 physicalEndIndex(QtyCentering centering, Direction direction) const;
    uint32 physicalEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                            Direction direction) const;
    template<typename NdArrayImpl>
    uint32 physicalEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                            Direction direction) const;

    uint32 ghostStartIndex(QtyCentering centering, Direction direction) const;
    uint32 ghostStartIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const;
    template<typename NdArrayImpl>
    uint32 ghostStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                           Direction direction) const;

    uint32 ghostEndIndex(QtyCentering centering, Direction direction) const;
    uint32 ghostEndIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const;
    template<typename NdArrayImpl>
    uint32 ghostEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                         Direction direction) const;

    std::array<uint32, dim> allocSize(HybridQuantity::Scalar qtyType) const;
    std::array<uint32, dim> allocSizeDerived(HybridQuantity::Scalar qty, Direction dir) const;

    template<typename NdArrayImpl>
    void deriv(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand, Direction direction,
               Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const;

    template<typename NdArrayImpl, typename... Indexes>
    Point<double, dim> fieldNodeCoordinates(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                                            const Point<double, dim>& origin,
                                            Indexes... index) const;

    template<typename... Indexes>
    Point<double, dim> cellCenteredCoordinates(Indexes... index) const;

    template<typename NdArrayImpl>
    QtyCentering fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                Direction dir) const;

    uint32 nbrGhostNodes(QtyCentering const& centering) const;
    template<typename NdArrayImpl>
    uint32 nbrGhostNodes(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                         Direction direction) const;

    std::array<uint32, dim> nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const;
    template<typename NdArrayImpl>
    std::array<uint32, dim>
    nbrPhysicalNodes(Field<NdArrayImpl, HybridQuantity::Scalar> const& field) const;



    // routines for projecting some quantities form one centering to another
    // not all quantities re-centering need to be coded, we only write
    // those needed for Ohm, namely those for vector product
    LinearCombination const& momentsToEx() const { return impl_.momentsToEx(); }
    LinearCombination const& momentsToEy() const { return impl_.momentsToEy(); }
    LinearCombination const& momentsToEz() const { return impl_.momentsToEz(); }

    LinearCombination const& ByToEx() const { return impl_.ByToEx(); }
    LinearCombination const& ByToEz() const { return impl_.ByToEz(); }

    LinearCombination const& BxToEy() const { return impl_.BxToEy(); }
    LinearCombination const& BxToEz() const { return impl_.BxToEz(); }

    LinearCombination const& BzToEx() const { return impl_.BzToEx(); }
    LinearCombination const& BzToEy() const { return impl_.BzToEy(); }

    LinearCombination const& ExToMoment() const { return impl_.ExToMoment(); }
    LinearCombination const& EyToMoment() const { return impl_.EyToMoment(); }
    LinearCombination const& EzToMoment() const { return impl_.EzToMoment(); }
};

/* ---------------------------------------------------------------------------
 *
 *                                  PUBLIC
 *
 * --------------------------------------------------------------------------- */

template<typename DataType>
constexpr DataType inverse(DataType const& data)
{
    return static_cast<DataType>(1.) / data;
}

/**
 * @brief GridLayout<layout,dim>::GridLayout constructs a GridLayout
 * @param dxdydz  mesh resolution in each direction. Used dimensions must
 * have a non-zero mesh size and vice-versa otherwise a runtime_error exception
 * is thrown.
 * @param nbrCells number of cells in the physical domain in each direction.
 * All used dimensions must have at least GridLayout<layout,dim>::minNbrCells=10 cells otherwise
 * a runtime_error exception is thrown. Each non-zero nbr Cell must be associated
 * with a non-zero mesh size (dxdydz) and vice-versa.
 * @param nbDims number of dimenions in the problem, can be 1, 2 or 3
 * @param layoutName is the name of the specific Grid Layout desired. For now only
 * 'yee' is available.
 * @param ghostParameter is an integer GridLayout uses to determine the number
 * of ghost nodes required for each quantity in each direction.
 */
template<Layout layout, std::size_t dim>
GridLayout<layout, dim>::GridLayout(std::array<double, dim> const& dx,
                                    std::array<uint32, dim> const& nbrCells,
                                    std::string const& layoutName, Point<double, dim> const& origin,
                                    uint32 ghostParameter)

    : interpOrder_{ghostParameter}
    , layoutName_{layoutName}
    , impl_{dx, nbrCells, origin, ghostParameter}
{
    origin_ = origin;

    dl_ = dx;
    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        odxyz_[iDir] = inverse(dl_[iDir]);
    }
    nbrCell_ = nbrCells;

    throwNotValid();
}




/**
 * @brief GridLayout<layout,dim>::physicalStartIndex return the index of the first primal or
 * dual node (centering= QtyCentering::dual or QtyCentering::primal).
 * The function just calls its private implementation.
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::physicalStartIndex(QtyCentering centering,
                                                   Direction direction) const
{
    return impl_.physicalStartIndex(centering, direction);
}


/**
 * @brief GridLayout<layout,dim>::physicalStartIndex return the index of the first node
 * of a 'HybridQuantity::Quantity' in a given 'direction' that belong to the physical domain.
 * The function just calls its private implementation
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::physicalStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                   Direction direction) const
{
    return impl_.physicalStartIndex(hybridQuantity, direction);
}


/**
 * @brief GridLayout::physicalStartIndex return the index of the first node
 * of a 'field' in a given 'direction' that belong to the physical domain.
 * The function just calls its private implementation
 */
template<Layout layout, std::size_t dim>
template<typename NdArrayImpl>
uint32
GridLayout<layout, dim>::physicalStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                            Direction direction) const
{
    return impl_.physicalStartIndex(field, direction);
}



/**
 * @brief GridLayout<layout,dim>::physicalEndIndex return the index of the last primal or
 * dual node in a given direction.
 * @param centering is either QtyCentering::dual or QtyCentering::primal
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::physicalEndIndex(QtyCentering centering, Direction direction) const
{
    return impl_.physicalEndIndex(centering, direction);
}


/**
 * @brief GridLayout<layout,dim>::physicalEndIndex return the index of the last node of a
 * 'HybridQuantity::Quantity' in a given 'direction'. The function just calls its private
 * implementation
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::physicalEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                 Direction direction) const
{
    return impl_.physicalEndIndex(hybridQuantity, direction);
}

/**
 * @brief GridLayout::physicalEndIndex return the index of the last node of a
 * 'field' in a given 'direction'. The function just calls its private implementation
 */
template<Layout layout, std::size_t dim>
template<typename NdArrayImpl>
uint32
GridLayout<layout, dim>::physicalEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                          Direction direction) const
{
    return impl_.physicalEndIndex(field, direction);
}



/**
 * @brief GridLayout<layout,dim>::ghostStartIndex returns the index of the first ghost node
 * of a 'field' in a given 'direction'
 * @param centering is either QtyCentering::dual or QtyCentering::primal
 * @param direction
 * @return
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::ghostStartIndex(QtyCentering centering, Direction direction) const
{
    return impl_.ghostStartIndex(centering, direction);
}

/**
 * @brief GridLayout<layout,dim>::ghostStartIndex returns the index of the first ghost node
 * of a 'HybridQuantity::Quantity' in a given 'direction'. The function just calls its private
 * implementation.
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::ghostStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                Direction direction) const
{
    return impl_.ghostStartIndex(hybridQuantity, direction);
}


/**
 * @brief GridLayout::ghostStartIndex returns the index of the first ghost node
 * of a 'field' in a given 'direction'. The function just calls its private
 * implementation.
 */
template<Layout layout, std::size_t dim>
template<typename NdArrayImpl>
uint32
GridLayout<layout, dim>::ghostStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                         Direction direction) const
{
    return impl_.ghostStartIndex(field, direction);
}

/**
 * @brief GridLayout<layout,dim>::ghostEndIndex returns the index of the last ghost node
 * of a 'field' in a given 'direction'
 * @param centering is either QtyCentering::dual or QtyCentering::primal
 * @param direction
 * @return
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::ghostEndIndex(QtyCentering centering, Direction direction) const
{
    return impl_.ghostEndIndex(centering, direction);
}

/**
 * @brief GridLayout<layout,dim>::ghostEndIndex returns the index of the last primal or dual
 * node in a given 'direction'.
 */
template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::ghostEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                                              Direction direction) const
{
    return impl_.ghostEndIndex(hybridQuantity, direction);
}


/**
 * @brief GridLayout::ghostEndIndex returns the index of the last primal or dual
 * node in a given 'direction'.
 */
template<Layout layout, std::size_t dim>
template<typename NdArrayImpl>
uint32
GridLayout<layout, dim>::ghostEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                       Direction direction) const
{
    return impl_.ghostEndIndex(field, direction);
}


/**
 * @brief GridLayout<layout,dim>::deriv calculate and return the derivative of a Field
 * in a given direction. In practice the function just calls its private implementation
 */
template<Layout layout, std::size_t dim>
template<typename NdArrayImpl>
void GridLayout<layout, dim>::deriv(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
                                    Direction direction,
                                    Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const
{
    switch (nbDims_)
    {
        case 1: impl_.deriv1D(operand, derivative); break;

        default: throw std::runtime_error("Error - only 1D fields are handled for the moment ");
    }
}




/**
 * @brief GridLayout<layout,dim>::allocSize
 * @return An std::array<uint32, dim> object, containing the allocation size of arrays in all
 * directions
 */
template<Layout layout, std::size_t dim>
std::array<uint32, dim> GridLayout<layout, dim>::allocSize(HybridQuantity::Scalar qtyType) const
{
    return impl_.allocSize(qtyType);
}




/**
 * @brief GridLayout<layout,dim>::allocSizeDerived calculate the allocation size in each
 * direction for a field that is the first derivative of an HybridQuantity::Quantity in
 * a given direction.
 * @param qty the quantity to be derivated
 * @param dir direciton of derivation
 * @return an std::array<uint32, dim> containing the size of the derivative array in each direction
 */
template<Layout layout, std::size_t dim>
std::array<uint32, dim> GridLayout<layout, dim>::allocSizeDerived(HybridQuantity::Scalar qty,
                                                                  Direction dir) const
{
    return impl_.allocSizeDerived(qty, dir);
}

/**
 * @brief GridLayout::fieldNodeCoordinates calculate the physical coordinate
 * of a certain mesh point (ix,iy,iz) for a given quantity.
 * @return a Point centered at the field node coordinate.
 */
template<Layout layout, std::size_t dim>
template<typename NdArrayImpl, typename... Indexes>
Point<double, dim> GridLayout<layout, dim>::fieldNodeCoordinates(
    const Field<NdArrayImpl, HybridQuantity::Scalar>& field, const Point<double, dim>& origin,
    Indexes... index) const
{
    return impl_.fieldNodeCoordinates(field, origin, index...);
}



/**
 * @brief GridLayout<layout,dim>::cellCenteredCoordinates calculates the physical coordinate
 * of a cell center indexed (ix,iy,iz).
 * @return a Point coordinate centered at the cell center.
 */
template<Layout layout, std::size_t dim>
template<typename... Indexes>
Point<double, dim> GridLayout<layout, dim>::cellCenteredCoordinates(Indexes... index) const
{
    return impl_.cellCenteredCoordinates(index...);
}

template<Layout layout, std::size_t dim>
template<typename NdArrayImpl>
QtyCentering
GridLayout<layout, dim>::fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                        Direction dir) const
{
    return impl_.fieldCentering(field, dir);
}



template<Layout layout, std::size_t dim>
uint32 GridLayout<layout, dim>::nbrGhostNodes(QtyCentering const& centering) const
{
    return impl_.nbrGhostNodes(centering);
}




template<Layout layout, std::size_t dim>
std::array<uint32, dim>
GridLayout<layout, dim>::nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const
{
    return impl_.nbrPhysicalNodes(hybQty);
}



/* ---------------------------------------------------------------------------
 *
 *                                  PRIVATE
 *
 * --------------------------------------------------------------------------- */


template<Layout layout, std::size_t dim>
const std::string GridLayout<layout, dim>::errorInverseMesh = "GridLayout error: Invalid use of";


template<Layout layout, std::size_t dim>
void GridLayout<layout, dim>::throwNotValid() const
{
    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        if (dl_[iDir] == 0.0)
            throw std::runtime_error("Error - requires non 0 for any component of dl");

        if (dl_[iDir] < 0)
            throw std::runtime_error("Error - dl must be positive");

        if (nbrCell_[iDir] < minNbrCells)
            throw std::runtime_error("Error - number cells is too small");
    }
}


} // namespace PHARE

#endif // PHARE_CORE_DATA_GRID_GRIDLAYOUT_H
