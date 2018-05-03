#ifndef PHARE_CORE_GRID_GRIDLAYOUTIMPLINTERNALS_H
#define PHARE_CORE_GRID_GRIDLAYOUTIMPLINTERNALS_H


#include <array>
#include <cmath>

#include "gridlayoutdefs.h"

#include "data/field/field.h"
#include "hybrid/hybrid_quantities.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"
#include "utilities/types.h"

namespace PHARE
{
/**
 * @brief GridLayoutImplInternals is intended to factorize attributes and methods
 * common to all GridLayoutImpl derived classes (ex: GridLayoutImplYee).
 *
 * Most of the implementations needed to handle GridLayout operations are provided
 * by GridLayoutImplInternals' methods. A lot of operations on a quantity
 * actually only depend on the centering of the quantity, namely whether it is
 * a primal or dual centering. The only thing GridlayoutImplInternals
 *
 */
template<typename GridImpl, std::size_t dim>
class GridLayoutImplInternals
{
    // ------------------------------------------------------------------------
    //                              PRIVATE
    //
    //   this is just for GridLayoutImplInternals
    // ------------------------------------------------------------------------
private:
    std::array<uint32, dim>
    nodeNbrFromCentering_(std::array<QtyCentering, dim> const& qtyCenterings) const;



    // ------------------------------------------------------------------------
    //                             PROTECTED
    //
    // this code will be shared by all concrete of GridLayoutImpl*
    // ------------------------------------------------------------------------
protected:
    /* uint32 static constexpr nbdims_{GridImpl::getDimensions()}; */
    uint32 constexpr static nbdims_{dim};

    uint32 nbrPrimalGhosts_;
    uint32 nbrDualGhosts_;
    uint32 interpOrder_;

    std::array<double, dim> dl_;

    Point<double, dim> origin_;

    std::array<double, dim> odxdydz_;
    std::array<uint32, dim> nbrPhysicalCells_;


    // stores key indices in each direction (3) for primal and dual nodes (2)
    std::array<std::array<uint32, 3>, 2> physicalStartIndexTable_;
    std::array<std::array<uint32, 3>, 2> physicalEndIndexTable_;
    std::array<std::array<uint32, 3>, 2> ghostEndIndexTable_;


    GridLayoutImplInternals(std::array<double, dim> const& dx,
                            std::array<uint32, dim> const& nbrCells, uint32 ghostParameter,
                            Point<double, dim> const& origin);




    double inverseSpatialStep(Direction direction) const noexcept
    {
        return odxdydz_[static_cast<uint32>(direction)];
    }


    void computeNbrGhosts(uint32 ghostParameter);

    std::array<uint32, dim>
    physicalNodeNbrFromCentering_(std::array<QtyCentering, dim> const& qtyCenterings) const;

    uint32 physicalStartIndex_(QtyCentering centering, Direction direction) const;
    uint32 physicalStartIndex_(HybridQuantity::Scalar const& hybridQuantity,
                               Direction direction) const;

    template<typename NdArrayImpl>
    uint32 physicalStartIndex_(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                               Direction direction) const;

    uint32 physicalEndIndex_(QtyCentering centering, Direction direction) const;
    uint32 physicalEndIndex_(HybridQuantity::Scalar const& hybridQuantity,
                             Direction direction) const;

    template<typename NdArrayImpl>
    uint32 physicalEndIndex_(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                             Direction direction) const;


    uint32 ghostStartIndex_(QtyCentering centering, Direction direction) const;
    uint32 ghostStartIndex_(HybridQuantity::Scalar const& hybridQuantity,
                            Direction direction) const;

    template<typename NdArrayImpl>
    uint32 ghostStartIndex_(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                            Direction direction) const;

    uint32 ghostEndIndex_(QtyCentering centering, Direction direction) const;
    uint32 ghostEndIndex_(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const;

    template<typename NdArrayImpl>
    uint32 ghostEndIndex_(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                          Direction direction) const;


    std::array<uint32, dim> allocSize_(HybridQuantity::Scalar qty) const;
    std::array<uint32, dim> allocSizeDerived_(HybridQuantity::Scalar qty, Direction dir) const;


    template<typename NdArrayImpl, typename... Indexes>
    Point<double, dim>
    fieldNodeCoordinates_(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                          const Point<double, dim>& origin, Indexes... index) const;

    template<typename... Indexes>
    Point<double, dim> cellCenteredCoordinates_(Indexes... index) const;

    template<typename NdArrayImpl>
    void deriv1D_(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
                  Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const;

    template<typename NdArrayImpl>
    QtyCentering fieldCentering_(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                 Direction dir) const;

    void initPhysicalStart(const gridDataT& data);
    void initPhysicalEnd(const gridDataT& data);
    void initGhostStart(const gridDataT& data);
    void initGhostEnd(const gridDataT& data);

    QtyCentering changeCentering(QtyCentering layout) const;
    QtyCentering derivedCentering(HybridQuantity::Scalar qty, Direction dir) const;


    uint32 nbrGhosts(QtyCentering centering) const noexcept;
    uint32 isDual(QtyCentering centering) const noexcept;

    uint32 ghostCellIndexAtMax(QtyCentering centering, Direction direction) const;

private:
    void initInternals_(std::array<double, dim> const& dxdydz,
                        std::array<uint32, dim> const& nbrCells, Point<double, dim> const& origin);
};
template<typename GridImpl, std::size_t dim>
GridLayoutImplInternals<GridImpl, dim>::GridLayoutImplInternals(
    std::array<double, dim> const& dx, std::array<uint32, dim> const& nbrCells,
    uint32 ghostParameter, Point<double, dim> const& origin)
    : interpOrder_{ghostParameter}
{
    computeNbrGhosts(ghostParameter);

    initInternals_(dx, nbrCells, origin);
}

template<typename GridImpl, std::size_t dim>
void GridLayoutImplInternals<GridImpl, dim>::initInternals_(std::array<double, dim> const& dxdydz,
                                                            std::array<uint32, dim> const& nbrCells,
                                                            Point<double, dim> const& origin)
{
    this->dl_     = dxdydz;
    this->origin_ = origin;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        odxdydz_[iDir] = 1. / dl_[iDir];
    }
    nbrPhysicalCells_ = nbrCells;
}
// TODO : mark gridDataT as dim dependent
template<typename GridImpl, std::size_t dim>
void GridLayoutImplInternals<GridImpl, dim>::initPhysicalStart(const gridDataT& data)
{
    uint32 iprimal = static_cast<uint32>(data.primal);
    uint32 idual   = static_cast<uint32>(data.dual);

    physicalStartIndexTable_[iprimal][data.idirX] = nbrGhosts(QtyCentering::primal);
    physicalStartIndexTable_[iprimal][data.idirY] = nbrGhosts(QtyCentering::primal);
    physicalStartIndexTable_[iprimal][data.idirZ] = nbrGhosts(QtyCentering::primal);

    physicalStartIndexTable_[idual][data.idirX] = nbrGhosts(QtyCentering::dual);
    physicalStartIndexTable_[idual][data.idirY] = nbrGhosts(QtyCentering::dual);
    physicalStartIndexTable_[idual][data.idirZ] = nbrGhosts(QtyCentering::dual);
}



/**
 * @brief GridLayoutImplInternals<GridImpl::dim>::initPhysicalEnd intialize the table of indices
 * corresponding to the last node for primal and dual centering.
 * The formula is simple : the last index is obtained from the first one
 * (which is physicalStartIndex of primal/dual in a given direction)
 *  + the number of cells minus 1 for dual nodes only.
 */
template<typename GridImpl, std::size_t dim>
void GridLayoutImplInternals<GridImpl, dim>::initPhysicalEnd(const gridDataT& data)
{
    uint32 iprimal = static_cast<uint32>(data.primal);
    uint32 idual   = static_cast<uint32>(data.dual);

    physicalEndIndexTable_[iprimal][data.idirX] = physicalStartIndexTable_[iprimal][data.idirX]
                                                  + nbrPhysicalCells_[data.idirX]
                                                  - isDual(data.primal);

    physicalEndIndexTable_[iprimal][data.idirY] = physicalStartIndexTable_[iprimal][data.idirY]
                                                  + nbrPhysicalCells_[data.idirY]
                                                  - isDual(data.primal);

    physicalEndIndexTable_[iprimal][data.idirZ] = physicalStartIndexTable_[iprimal][data.idirZ]
                                                  + nbrPhysicalCells_[data.idirZ]
                                                  - isDual(data.primal);


    physicalEndIndexTable_[idual][data.idirX] = physicalStartIndexTable_[idual][data.idirX]
                                                + nbrPhysicalCells_[data.idirX] - isDual(data.dual);

    physicalEndIndexTable_[idual][data.idirY] = physicalStartIndexTable_[idual][data.idirY]
                                                + nbrPhysicalCells_[data.idirY] - isDual(data.dual);

    physicalEndIndexTable_[idual][data.idirZ] = physicalStartIndexTable_[idual][data.idirZ]
                                                + nbrPhysicalCells_[data.idirZ] - isDual(data.dual);
}




/**
 * @brief GridLayoutImplInternals<GridImpl::dim>::initGhostEnd calculate and stores the index
 * of the last primal and dual nodes in each direction. The formula simply
 * consists in starting at physicalEndIndex() and to add the number of ghost nodes.
 */
template<typename GridImpl, std::size_t dim>
void GridLayoutImplInternals<GridImpl, dim>::initGhostEnd(const gridDataT& data)
{
    uint32 iprimal = static_cast<uint32>(data.primal);
    uint32 idual   = static_cast<uint32>(data.dual);

    ghostEndIndexTable_[iprimal][data.idirX]
        = physicalEndIndexTable_[iprimal][data.idirX] + nbrGhosts(data.primal);

    ghostEndIndexTable_[iprimal][data.idirY]
        = physicalEndIndexTable_[iprimal][data.idirY] + nbrGhosts(data.primal);

    ghostEndIndexTable_[iprimal][data.idirZ]
        = physicalEndIndexTable_[iprimal][data.idirZ] + nbrGhosts(data.primal);

    ghostEndIndexTable_[idual][data.idirX]
        = physicalEndIndexTable_[idual][data.idirX] + nbrGhosts(data.dual);

    ghostEndIndexTable_[idual][data.idirY]
        = physicalEndIndexTable_[idual][data.idirY] + nbrGhosts(data.dual);

    ghostEndIndexTable_[idual][data.idirZ]
        = physicalEndIndexTable_[idual][data.idirZ] + nbrGhosts(data.dual);
}




/**
 * @brief GridLayoutImplInternals<GridImpl::dim>::derivedCentering this function returns the
 * centering (primal or dual) of a quantity after a first order derivation. dual becomes primal and
 * primal becomes dual. hybridQuantityCentering is used to know if the HybridQuantity::Quantity
 * 'qty' is primal or dual in the Direction 'dir'
 */
template<typename GridImpl, std::size_t dim>
QtyCentering GridLayoutImplInternals<GridImpl, dim>::derivedCentering(HybridQuantity::Scalar qty,
                                                                      Direction dir) const
{
    uint32 iField = static_cast<uint32>(qty);
    uint32 idir   = static_cast<uint32>(dir);


    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;

    QtyCentering newCentering = changeCentering(hybridQtyCentering[iField][idir]);

    return newCentering;
}


template<typename GridImpl, std::size_t dim>
std::array<uint32, dim> GridLayoutImplInternals<GridImpl, dim>::physicalNodeNbrFromCentering_(
    std::array<QtyCentering, dim> const& qtyCenterings) const
{
    std::array<uint32, dim> nodeNbr;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        nodeNbr[iDir] = nbrPhysicalCells_[iDir] + 1 - isDual(qtyCenterings[iDir]);
    }

    return nodeNbr;
}




/**
 * @brief GridLayoutImplInternals<GridImpl::dim>::nodeNbrFromCentering_ returns an array containing
 * the total number of nodes (ghosts + physical) in each direction.
 * The calculation is easy : there are nbrPhysicalCells + 1 nodes in the domain
 * + 2 times the number of ghost nodes.
 */
template<typename GridImpl, std::size_t dim>
std::array<uint32, dim> GridLayoutImplInternals<GridImpl, dim>::nodeNbrFromCentering_(
    std::array<QtyCentering, dim> const& qtyCenterings) const
{
    std::array<uint32, dim> nbrNodes = physicalNodeNbrFromCentering_(qtyCenterings);

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        nbrNodes[iDir] += 2 * nbrGhosts(qtyCenterings[iDir]);
    }


    return nbrNodes;
}




/**
 * @brief GridLayoutImplInternals<GridImpl::dim>::allocSize_
 * @return An std::array<uint32, dim> object, containing the size to which allocate arrays
 * of an HybridQuantity::Quantity 'qty' in every directions.
 */
template<typename GridImpl, std::size_t dim>
std::array<uint32, dim>
GridLayoutImplInternals<GridImpl, dim>::allocSize_(HybridQuantity::Scalar qty) const
{
    uint32 iQty = static_cast<uint32>(qty);


    // TODO: hybridQtyCentering should be defined per dimension so that we could simply do
    // auto sizeArray = nodeNbrFromCentering_(hybridQtyCentering[iQty]);

    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;

    std::array<QtyCentering, dim> qtyCentering;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        qtyCentering[iDir] = hybridQtyCentering[iQty][iDir];
    }


    return nodeNbrFromCentering_(qtyCentering);
}




template<typename GridImpl, std::size_t dim>
std::array<uint32, dim>
GridLayoutImplInternals<GridImpl, dim>::allocSizeDerived_(HybridQuantity::Scalar qty,
                                                          Direction dir) const
{
    uint32 iDerivedDir = static_cast<uint32>(dir);
    uint32 iQty        = static_cast<uint32>(qty);

    // get the centering of the derivative of 'qty' in the direction of derivation
    QtyCentering newCentering = derivedCentering(qty, dir);



    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;

    std::array<QtyCentering, dim> qtyCenterings;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        qtyCenterings[iDir] = hybridQtyCentering[iQty][iDir];
    }




    // ...and permute the centering in the direction of derivation
    qtyCenterings[iDerivedDir] = newCentering;

    // get the total number of nodes (ghost + physical) for the new centering
    return nodeNbrFromCentering_(qtyCenterings);
}

template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::physicalStartIndex_(QtyCentering centering,
                                                                   Direction direction) const
{
    uint32 icentering = static_cast<uint32>(centering);
    uint32 iDir       = static_cast<uint32>(direction);
    return physicalStartIndexTable_[icentering][iDir];
}


template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::physicalStartIndex_(
    HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
{
    uint32 iQty                       = static_cast<uint32>(hybridQuantity);
    uint32 iDir                       = static_cast<uint32>(direction);
    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;
    uint32 iCentering                 = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);

    return physicalStartIndexTable_[iCentering][iDir];
}
template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImplInternals<GridImpl, dim>::physicalStartIndex_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    return physicalStartIndex_(field.physicalQuantity(), direction);
}



template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::physicalEndIndex_(QtyCentering centering,
                                                                 Direction direction) const
{
    uint32 icentering = static_cast<uint32>(centering);
    uint32 iDir       = static_cast<uint32>(direction);

    return physicalEndIndexTable_[icentering][iDir];
}

template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::physicalEndIndex_(
    HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
{
    uint32 iQty                       = static_cast<uint32>(hybridQuantity);
    uint32 iDir                       = static_cast<uint32>(direction);
    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;
    uint32 iCentering                 = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);

    return physicalEndIndexTable_[iCentering][iDir];
}

template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImplInternals<GridImpl, dim>::physicalEndIndex_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    return physicalEndIndex_(field.physicalQuantity(), direction);
}


template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::ghostStartIndex_(QtyCentering centering,
                                                                Direction direction) const
{
    // ghostStartIndex is always the first node
    return 0;
}


template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::ghostStartIndex_(
    HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
{
    // ghostStartIndex is always the first node
    return 0;
}


template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImplInternals<GridImpl, dim>::ghostStartIndex_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    // ghostStartIndex is always the first node
    return 0;
}




template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::ghostEndIndex_(QtyCentering centering,
                                                              Direction direction) const
{
    uint32 iCentering = static_cast<uint32>(centering);
    uint32 iDir       = static_cast<uint32>(direction);

    return ghostEndIndexTable_[iCentering][iDir];
}

template<typename GridImpl, std::size_t dim>
uint32
GridLayoutImplInternals<GridImpl, dim>::ghostEndIndex_(HybridQuantity::Scalar const& hybridQuantity,
                                                       Direction direction) const
{
    uint32 iQty                       = static_cast<uint32>(hybridQuantity);
    uint32 iDir                       = static_cast<uint32>(direction);
    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;
    uint32 iCentering                 = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);
    return ghostEndIndexTable_[iCentering][iDir];
}

template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImplInternals<GridImpl, dim>::ghostEndIndex_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    return ghostEndIndex_(field.physicalQuantity(), direction);
}

template<typename GridImpl, std::size_t dim>
QtyCentering GridLayoutImplInternals<GridImpl, dim>::changeCentering(QtyCentering centering) const
{
    QtyCentering newCentering = QtyCentering::primal;

    if (centering == QtyCentering::primal)
    {
        newCentering = QtyCentering::dual;
    }

    return newCentering;
}



/**
 * @brief fieldNodeCoordinates returns a Point,
 * the idea is to make in every initializer
 * method, 3 nested loops over primal PhysicalStart/End indices and to get
 * the physical coordinate of those mesh nodes for the considered field.
 * @param field
 * the returned point depends on the field's centering
 * @param origin
 * @param ix is a primal or dual index
 * @param iy is a primal or dual index
 * @param iz is a primal or dual index
 * @return Point
 * the desired field-centered coordinate
 */
template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl, typename... Indexes>
Point<double, dim> GridLayoutImplInternals<GridImpl, dim>::fieldNodeCoordinates_(
    const Field<NdArrayImpl, HybridQuantity::Scalar>& field, const Point<double, dim>& origin,
    Indexes... index) const
{
    static_assert(sizeof...(Indexes) == dim, "Error dimension does not match number of arguments");


    uint32 iQuantity       = static_cast<uint32>(field.physicalQuantity());
    constexpr uint32 iDual = static_cast<uint32>(QtyCentering::dual);


    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;

    Point<int32, dim> coord{static_cast<int32>(index)...};

    Point<double, dim> position;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        double halfCell = 0.0;

        auto const centering = static_cast<uint32>(hybridQtyCentering[iQuantity][iDir]);
        int32 const iStart   = physicalStartIndexTable_[centering][iDir];

        // A shift of +dx/2, +dy/2, +dz/2 is necessary to get the physical
        // coordinate on the dual mesh
        // No shift for coordinate on the primal mesh
        // This shift DOES NOT DEPEND ON the interpolation order
        // Because,
        // if ix is primal then ixStart is primal
        // if ix is dual   then ixStart is dual
        // if iy is primal then iyStart is primal ...

        if (centering == iDual)
        {
            halfCell = 0.5;
        }


        position[iDir]
            = (static_cast<double>(coord[iDir] - iStart) + halfCell) * dl_[iDir] + origin[iDir];
    }

    return position;
}




/**
 * @brief cellCenteredCoordinates returns a cell-centered Point.
 * The idea is to call this method in every initializer method
 * using 3 nested loops over primal PhysicalStart/End indices.
 * This function will typically used to evaluate the density at initialization
 * phases, since it is usually defined at cell centers.
 * @param origin
 * @param ix is a primal index
 * @param iy is a primal index
 * @param iz is a primal index
 * @return Point
 * the desired cell-centered (dual/dual/dual) coordinate
 */
template<typename GridImpl, std::size_t dim>
template<typename... Indexes>
Point<double, dim>
GridLayoutImplInternals<GridImpl, dim>::cellCenteredCoordinates_(Indexes... index) const
{
    static_assert(sizeof...(Indexes) == dim, "Error dimension does not match number of arguments");

    uint32 constexpr iPrimal = static_cast<uint32>(QtyCentering::primal);

    constexpr double halfCell = 0.5;
    // A shift of +dx/2, +dy/2, +dz/2 is necessary to get the
    // cell center physical coordinates,
    // because this point is located on the dual mesh

    Point<uint32, dim> coord(index...);

    Point<double, dim> physicalPosition;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        auto iStart = physicalStartIndexTable_[iPrimal][iDir];

        physicalPosition[iDir]
            = (static_cast<double>(coord[iDir] - iStart) + halfCell) * dl_[iDir] + origin_[iDir];
    }

    return physicalPosition;




    // uint32 idirX   = static_cast<uint32>(Direction::X);
    // uint32 idirY   = static_cast<uint32>(Direction::Y);
    // uint32 idirZ   = static_cast<uint32>(Direction::Z);
    // uint32 iprimal = static_cast<uint32>(QtyCentering::primal);

    // uint32 ixStart = physicalStartIndexTable_[iprimal][idirX];
    // uint32 iyStart = physicalStartIndexTable_[iprimal][idirY];
    // uint32 izStart = physicalStartIndexTable_[iprimal][idirZ];


    // double x = ((ix - ixStart) + halfCell) * dx_ + origin_.x;
    // double y = ((iy - iyStart) + halfCell) * dy_ + origin_.y;
    // double z = ((iz - izStart) + halfCell) * dz_ + origin_.z;

    // return Point<double, dim>(physicalPosition);
}



template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl>
QtyCentering GridLayoutImplInternals<GridImpl, dim>::fieldCentering_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction dir) const
{
    uint32 iDir = static_cast<uint32>(dir);

    uint32 iQty                       = static_cast<uint32>(field.physicalQuantity());
    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;

    return hybridQtyCentering[iQty][iDir];
}


/**
 * @brief computeOffsets computes the number of ghost cells for fields.
 * On the primal mesh the number of ghosts depends on centeredOffset_
 * On the dual mesh the number of ghosts depends on
 * leftOffset_ and rightOffset_
 *
 * This is explained in details on the wiki page
 * https://hephaistos.lpp.polytechnique.fr/redmine/projects/hyb-par/wiki/PrimalDual
 *
 * @param ghostParameter, corresponds to the interpolation order
 */
template<typename GridImpl, std::size_t dim>
void GridLayoutImplInternals<GridImpl, dim>::computeNbrGhosts(uint32 ghostParameter)
{
    /* for first order Interpolation, there is no primal ghost node neeeded
       for particle/mesh interactions. However one ghost node is required
       for calculating Laplacians so we add one.
    */
    if (ghostParameter == 1)
    {
        nbrPrimalGhosts_ = static_cast<uint32>(std::floor((ghostParameter + 1) / 2.));
        nbrDualGhosts_   = static_cast<uint32>(std::floor((ghostParameter + 1) / 2.));
    }

    /* for interpolation order larger than 1, there is at least 1 primal ghost
       node so Laplacians can be calculated OK */
    else if (ghostParameter > 1)
    {
        nbrPrimalGhosts_ = static_cast<uint32>(std::floor(ghostParameter / 2.));
        nbrDualGhosts_   = static_cast<uint32>(std::floor((ghostParameter + 1) / 2.));
    }
}




/**
 * @brief GridLayoutImplInternals<GridImpl, dim>::nbrGhosts
 * @param centering QtyCentering::primal or QtyCentering::dual
 * @return the number of ghost nodes on each side of the mesh for a given centering
 */
template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::nbrGhosts(QtyCentering centering) const noexcept
{
    uint32 nbrGhosts = nbrPrimalGhosts_;

    if (centering == QtyCentering::dual)
    {
        nbrGhosts = nbrDualGhosts_;
    }
    return nbrGhosts;
}




/**
 * @brief GridLayoutImplInternals<GridImpl, dim>::isDual the method is used in index calculations
 * when we sometimes need to substract 1 when the centering is dual
 * @param centering QtyCentering::primal or QtyCentering::dual
 * @return return 1 if the centering is dual or 0 otherwise.
 */
template<typename GridImpl, std::size_t dim>
uint32 GridLayoutImplInternals<GridImpl, dim>::isDual(QtyCentering centering) const noexcept
{
    uint32 isdual = 0;

    if (centering == QtyCentering::dual)
    {
        isdual = 1;
    }

    return isdual;
}


/**
 * @brief GridLayoutImplYee::deriv1D It was decided to compute the
 * derivative on the entire physical domain.
 * In the case of Maxwell Ampere (dual centering of the operand),
 * it will therefore be necessary to get the values
 * of the operand outside the physical domain before applying
 * Maxwell Ampere
 *
 * @param operand is always primal in the case of Maxwell Faraday
 * rotational operator (dEz/dy, dEy/dz, dEx/dz, dEz/dx, dEy/dx, dEx/dy)
 * operand is always dual in the case of Maxwell Ampere
 * rotational operator (dBz/dy, dBy/dz, dBx/dz, dBz/dx, dBy/dx, dBx/dy)
 *
 * @param derivative
 */
template<typename GridImpl, std::size_t dim>
template<typename NdArrayImpl>
void GridLayoutImplInternals<GridImpl, dim>::deriv1D_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
    Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const
{
    uint32 iDirX = static_cast<uint32>(Direction::X);

    uint32 iQtyOperand                = static_cast<uint32>(operand.physicalQuantity());
    auto constexpr hybridQtyCentering = GridImpl::hybridQtyCentering_;

    QtyCentering opCentering = hybridQtyCentering[iQtyOperand][iDirX];


    // The QtyCentering of derivative is given by
    // iQty = static_cast<uint32>( derivative.physicalQuantity() )
    // hybridQtyCentering_[iQty][idir]
    uint32 iDerStart = physicalStartIndex_(derivative, Direction::X);
    uint32 iDerEnd   = physicalEndIndex_(derivative, Direction::X);

    uint32 iOpStart = physicalStartIndex_(operand, Direction::X);

    uint32 iOp = iOpStart;
    if (opCentering == QtyCentering::dual)
    {
        --iOp;
    }

    for (uint32 iDer = iDerStart; iDer <= iDerEnd; ++iDer)
    {
        derivative(iDer) = odxdydz_[0] * (operand(iOp + 1) - operand(iOp));
        ++iOp;
    }
}
} // namespace PHARE

#endif // GRIDLAYOUTIMPLINTERNALS_H
