#ifndef PHARE_CORE_GRID_GridLayout_H
#define PHARE_CORE_GRID_GridLayout_H


#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "data/field/field.h"
#include "gridlayoutdefs.h"
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
 * @brief GridLayout is intended to factorize attributes and methods
 * common to all GridLayoutImpl derived classes (ex: GridLayoutImplYee).
 *
 * Most of the implementations needed to handle GridLayout operations are provided
 * by GridLayout' methods. A lot of operations on a quantity
 * actually only depend on the centering of the quantity, namely whether it is
 * a primal or dual centering. The only thing GridLayout
 *
 */
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
class GridLayout
{
    // ------------------------------------------------------------------------
    //                             PROTECTED
    //
    // this code will be shared by all concrete of GridLayoutImpl*
    // ------------------------------------------------------------------------
public:
    /* uint32 static constexpr nbdims_{GridLayoutImpl::getDimensions()}; */
    // uint32 constexpr static nbdims_{dim};
    std::array<double, dim> meshSize_;
    Point<double, dim> origin_;
    std::array<uint32, dim> nbrPhysicalCells_;

    static constexpr gridDataT data{};



    // stores key indices in each direction (3) for primal and dual nodes (2)
    std::array<std::array<uint32, dim>, 2> physicalStartIndexTable_;
    std::array<std::array<uint32, dim>, 2> physicalEndIndexTable_;
    std::array<std::array<uint32, dim>, 2> ghostEndIndexTable_;

    std::array<double, dim> odxdydz_;

    GridLayout(std::array<double, dim> const& meshSize, std::array<uint32, dim> const& nbrCells,
               Point<double, dim> const& origin)
        : meshSize_{meshSize}
        , origin_{origin}
        , nbrPhysicalCells_{nbrCells}
        , physicalStartIndexTable_{initPhysicalStart()}
        , physicalEndIndexTable_{initPhysicalEnd()}
        , ghostEndIndexTable_{initGhostEnd()}
    {
    }

    auto initPhysicalStart()
    {
        std::array<std::array<uint32, dim>, 2> physicalStartIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        physicalStartIndexTable[iprimal][data.idirX] = nbrPrimalGhosts();
        physicalStartIndexTable[idual][data.idirX]   = nbrDualGhosts();

        if constexpr (dim > 1)
        {
            physicalStartIndexTable[iprimal][data.idirY] = nbrPrimalGhosts();
            physicalStartIndexTable[idual][data.idirY]   = nbrDualGhosts();

            if constexpr (dim > 2)
            {
                physicalStartIndexTable[iprimal][data.idirZ] = nbrPrimalGhosts();
                physicalStartIndexTable[idual][data.idirZ]   = nbrDualGhosts();
            }
        }
        return physicalStartIndexTable;
    }


    /**
     * @brief GridLayout<GridLayoutImpl::dim>::initPhysicalEnd intialize the table of indices
     * corresponding to the last node for primal and dual centering.
     * The formula is simple : the last index is obtained from the first one
     * (which is physicalStartIndex of primal/dual in a given direction)
     *  + the number of cells minus 1 for dual nodes only.
     */
    auto initPhysicalEnd()
    {
        std::array<std::array<uint32, dim>, 2> physicalEndIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        physicalEndIndexTable[iprimal][data.idirX]
            = physicalStartIndexTable_[iprimal][data.idirX] + nbrPhysicalCells_[data.idirX];


        physicalEndIndexTable[idual][data.idirX] = physicalStartIndexTable_[idual][data.idirX]
                                                   + nbrPhysicalCells_[data.idirX] - dualOffset();

        if constexpr (dim > 1)
        {
            physicalEndIndexTable[iprimal][data.idirY]
                = physicalStartIndexTable_[iprimal][data.idirY] + nbrPhysicalCells_[data.idirY];

            physicalEndIndexTable[idual][data.idirY] = physicalStartIndexTable_[idual][data.idirY]
                                                       + nbrPhysicalCells_[data.idirY]
                                                       - dualOffset();

            if constexpr (dim > 2)
            {
                physicalEndIndexTable[iprimal][data.idirZ]
                    = physicalStartIndexTable_[iprimal][data.idirZ] + nbrPhysicalCells_[data.idirZ];

                physicalEndIndexTable[idual][data.idirZ]
                    = physicalStartIndexTable_[idual][data.idirZ] + nbrPhysicalCells_[data.idirZ]
                      - dualOffset();
            }
        }
        return physicalEndIndexTable;
    }




    /**
     * @brief GridLayout<GridLayoutImpl::dim>::initGhostEnd calculate and stores the index
     * of the last primal and dual nodes in each direction. The formula simply
     * consists in starting at physicalEndIndex() and to add the number of ghost nodes.
     */
    auto initGhostEnd()
    {
        std::array<std::array<uint32, dim>, 2> ghostEndIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        ghostEndIndexTable[iprimal][data.idirX]
            = physicalEndIndexTable_[iprimal][data.idirX] + nbrPrimalGhosts();

        ghostEndIndexTable[idual][data.idirX]
            = physicalEndIndexTable_[idual][data.idirX] + nbrDualGhosts();

        if constexpr (dim > 1)
        {
            ghostEndIndexTable[iprimal][data.idirY]
                = physicalEndIndexTable_[iprimal][data.idirY] + nbrPrimalGhosts();

            ghostEndIndexTable[idual][data.idirY]
                = physicalEndIndexTable_[idual][data.idirY] + nbrDualGhosts();

            if constexpr (dim > 2)
            {
                ghostEndIndexTable[iprimal][data.idirZ]
                    = physicalEndIndexTable_[iprimal][data.idirZ] + nbrPrimalGhosts();

                ghostEndIndexTable[idual][data.idirZ]
                    = physicalEndIndexTable_[idual][data.idirZ] + nbrDualGhosts();
            }
        }
        return ghostEndIndexTable;
    }




    /**
     * @brief origin return the lower point of the grid described by the GridLayout
     * in physical coordinates
     */
    Point<double, dim> origin() const { return origin_; }


    /**
     * @brief returns the mesh size in the 'dim' dimensions
     */
    std::array<double, dim> meshSize() const { return meshSize_; }


    double inverseMeshSize(Direction direction) const noexcept
    {
        return odxdydz_[static_cast<uint32>(direction)];
    }


    std::array<double, dim> inverseMeshSize() const noexcept { return odxdydz_; }


    std::array<uint32, dim> nbrCells() const { return nbrPhysicalCells_; }



    uint32 nbDimensions() const { return dim; }



    // uint32 order() const { return interpOrder; }




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
    template<std::size_t order                         = interpOrder,
             std::enable_if_t<order == 1, dummy::type> = dummy::value>
    uint32 constexpr nbrDualGhosts() const
    {
        /* for first order Interpolation, there is no primal ghost node neeeded
           for particle/mesh interactions. However one ghost node is required
           for calculating Laplacians so we add one.
        */
        // nbrPrimalGhosts_ = static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
        // nbrDualGhosts_   = static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
        return static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
    }


    template<std::size_t order                         = interpOrder,
             std::enable_if_t<order == 1, dummy::type> = dummy::value>
    uint32 constexpr nbrPrimalGhosts() const
    {
        /* for first order Interpolation, there is no primal ghost node neeeded
           for particle/mesh interactions. However one ghost node is required
           for calculating Laplacians so we add one.
        */
        // nbrPrimalGhosts_ = static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
        // nbrDualGhosts_   = static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
        return static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
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
    template<std::size_t order                          = interpOrder,
             std::enable_if_t<(order > 1), dummy::type> = dummy::value>
    uint32 constexpr nbrDualGhosts() const
    {
        /* for interpolation order larger than 1, there is at least 1 primal ghost
           node so Laplacians can be calculated OK */
        return static_cast<uint32>(std::floor((interpOrder + 1) / 2.));
    }

    template<std::size_t order                          = interpOrder,
             std::enable_if_t<(order > 1), dummy::type> = dummy::value>
    uint32 constexpr nbrPrimalGhosts() const
    {
        /* for interpolation order larger than 1, there is at least 1 primal ghost
           node so Laplacians can be calculated OK */
        return static_cast<uint32>(std::floor(interpOrder / 2.));
        ;
    }




    uint32 physicalStartIndex(QtyCentering centering, Direction direction) const
    {
        uint32 icentering = static_cast<uint32>(centering);
        uint32 iDir       = static_cast<uint32>(direction);
        return physicalStartIndexTable_[icentering][iDir];
    }


    uint32 physicalStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                              Direction direction) const
    {
        uint32 iQty                       = static_cast<uint32>(hybridQuantity);
        uint32 iDir                       = static_cast<uint32>(direction);
        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
        uint32 iCentering                 = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);

        return physicalStartIndexTable_[iCentering][iDir];
    }


    template<typename NdArrayImpl>
    uint32 physicalStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                              Direction direction) const
    {
        return physicalStartIndex(field.physicalQuantity(), direction);
    }




    uint32 physicalEndIndex(QtyCentering centering, Direction direction) const
    {
        uint32 icentering = static_cast<uint32>(centering);
        uint32 iDir       = static_cast<uint32>(direction);

        return physicalEndIndexTable_[icentering][iDir];
    }



    uint32 physicalEndIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
    {
        uint32 iQty                       = static_cast<uint32>(hybridQuantity);
        uint32 iDir                       = static_cast<uint32>(direction);
        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
        uint32 iCentering                 = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);

        return physicalEndIndexTable_[iCentering][iDir];
    }



    template<typename NdArrayImpl>
    uint32 physicalEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                            Direction direction) const
    {
        return physicalEndIndex(field.physicalQuantity(), direction);
    }




    uint32 ghostStartIndex(QtyCentering centering, Direction direction) const
    {
        // ghostStartIndex is always the first node
        return 0;
    }



    uint32 ghostStartIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
    {
        // ghostStartIndex is always the first node
        return 0;
    }



    template<typename NdArrayImpl>
    uint32 ghostStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                           Direction direction) const
    {
        // ghostStartIndex is always the first node
        return 0;
    }



    uint32 ghostEndIndex(QtyCentering centering, Direction direction) const
    {
        uint32 iCentering = static_cast<uint32>(centering);
        uint32 iDir       = static_cast<uint32>(direction);

        return ghostEndIndexTable_[iCentering][iDir];
    }




    uint32 ghostEndIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
    {
        uint32 iQty                       = static_cast<uint32>(hybridQuantity);
        uint32 iDir                       = static_cast<uint32>(direction);
        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
        uint32 iCentering                 = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);
        return ghostEndIndexTable_[iCentering][iDir];
    }




    template<typename NdArrayImpl>
    uint32 ghostEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                         Direction direction) const
    {
        return ghostEndIndex(field.physicalQuantity(), direction);
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
    template<typename NdArrayImpl, typename... Indexes>
    Point<double, dim> fieldNodeCoordinates(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                                            const Point<double, dim>& origin,
                                            Indexes... index) const
    {
        static_assert(sizeof...(Indexes) == dim,
                      "Error dimension does not match number of arguments");


        uint32 iQuantity       = static_cast<uint32>(field.physicalQuantity());
        constexpr uint32 iDual = static_cast<uint32>(QtyCentering::dual);


        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

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
                = (static_cast<double>(coord[iDir] - iStart) + halfCell) * meshSize_[iDir]
                  + origin[iDir];
        }

        return position;
    }




    template<typename... Indexes>
    Point<double, dim> cellCenteredCoordinates(Indexes... index) const
    {
        static_assert(sizeof...(Indexes) == dim,
                      "Error dimension does not match number of arguments");

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
                = (static_cast<double>(coord[iDir] - iStart) + halfCell) * meshSize_[iDir]
                  + origin_[iDir];
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


    /**
     * @brief GridLayout<GridLayoutImpl, dim>::nbrGhosts
     * @param centering QtyCentering::primal or QtyCentering::dual
     * @return the number of ghost nodes on each side of the mesh for a given centering
     */
    uint32 constexpr nbrGhosts(QtyCentering centering) const
    {
        uint32 nbrGhosts = nbrPrimalGhosts();

        if (centering == QtyCentering::dual)
        {
            nbrGhosts = nbrDualGhosts();
        }
        return nbrGhosts;
    }




    QtyCentering changeCentering(QtyCentering centering) const
    {
        QtyCentering newCentering = QtyCentering::primal;

        if (centering == QtyCentering::primal)
        {
            newCentering = QtyCentering::dual;
        }

        return newCentering;
    }



    uint32 constexpr dualOffset() const noexcept { return 1; }




    // ----------------------------------------------------------------------
    //                      LAYOUT SPECIFIC METHODS
    //
    // The methods below return results that need information specific to the
    // layout that is used
    // ----------------------------------------------------------------------


    std::string layoutName() const { return GridLayoutImpl::layoutName_; }


    /**
     * @brief returns the centering of all hybrid quantities in all directions
     */
    constexpr static std::array<QtyCentering, dim>
    centering(HybridQuantity::Scalar const& hybridQuantity)
    {
        return GridLayoutImpl::centering(hybridQuantity);
    }




    /**
     * @brief GridLayout<GridLayoutImpl::dim>::allocSize_
     * @return An std::array<uint32, dim> object, containing the size to which allocate arrays
     * of an HybridQuantity::Quantity 'qty' in every directions.
     */
    std::array<uint32, dim> allocSize(HybridQuantity::Scalar qty) const
    {
        uint32 iQty = static_cast<uint32>(qty);


        // TODO: hybridQtyCentering should be defined per dimension so that we could simply do
        // auto sizeArray = nodeNbrFromCentering_(hybridQtyCentering[iQty]);

        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        std::array<QtyCentering, dim> qtyCentering;

        for (std::size_t iDir = 0; iDir < dim; ++iDir)
        {
            qtyCentering[iDir] = hybridQtyCentering[iQty][iDir];
        }

        return nodeNbrFromCentering_(qtyCentering);
    }



    std::array<uint32, dim> allocSizeDerived(HybridQuantity::Scalar qty, Direction dir) const
    {
        uint32 iDerivedDir = static_cast<uint32>(dir);
        uint32 iQty        = static_cast<uint32>(qty);

        // get the centering of the derivative of 'qty' in the direction of derivation
        QtyCentering newCentering = derivedCentering(qty, dir);

        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

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




    template<typename NdArrayImpl>
    void deriv1D_(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
                  Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const;

    template<typename NdArrayImpl>
    QtyCentering fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                Direction dir) const
    {
        uint32 iDir                       = static_cast<uint32>(dir);
        uint32 iQty                       = static_cast<uint32>(field.physicalQuantity());
        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        return hybridQtyCentering[iQty][iDir];
    }




    std::array<uint32, dim> nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const
    {
        std::array<QtyCentering, dim> centerings;

        for (std::size_t iDir = 0; iDir < dim; ++iDir)
        {
            centerings[iDir]
                = GridLayoutImpl::hybridQtyCentering_[static_cast<uint32>(hybQty)][iDir];
        }

        return this->physicalNodeNbrFromCentering_(centerings);
    }



    /**
     * @brief GridLayout<GridLayoutImpl::dim>::derivedCentering this function returns the
     * centering (primal or dual) of a quantity after a first order derivation. dual becomes primal
     * and primal becomes dual. hybridQuantityCentering is used to know if the
     * HybridQuantity::Quantity 'qty' is primal or dual in the Direction 'dir'
     */
    QtyCentering derivedCentering(HybridQuantity::Scalar qty, Direction dir) const
    {
        uint32 iField = static_cast<uint32>(qty);
        uint32 idir   = static_cast<uint32>(dir);


        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        QtyCentering newCentering = changeCentering(hybridQtyCentering[iField][idir]);

        return newCentering;
    }




private:
    std::array<uint32, dim>
    physicalNodeNbrFromCentering_(std::array<QtyCentering, dim> const& qtyCenterings) const
    {
        std::array<uint32, dim> nodeNbr;

        for (std::size_t iDir = 0; iDir < dim; ++iDir)
        {
            nodeNbr[iDir] = nbrPhysicalCells_[iDir] + 1
                            - ((qtyCenterings[iDir] == QtyCentering::dual) ? dualOffset() : 0);
        }

        return nodeNbr;
    }




    /**
     * @brief GridLayout<GridLayoutImpl::dim>::nodeNbrFromCentering_ returns an array containing
     * the total number of nodes (ghosts + physical) in each direction.
     * The calculation is easy : there are nbrPhysicalCells + 1 nodes in the domain
     * + 2 times the number of ghost nodes.
     */
    std::array<uint32, dim>
    nodeNbrFromCentering_(std::array<QtyCentering, dim> const& qtyCenterings) const
    {
        std::array<uint32, dim> nbrNodes = physicalNodeNbrFromCentering_(qtyCenterings);

        for (std::size_t iDir = 0; iDir < dim; ++iDir)
        {
            nbrNodes[iDir] += 2 * nbrGhosts(qtyCenterings[iDir]);
        }

        return nbrNodes;
    }
};




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
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
template<typename NdArrayImpl>
void GridLayout<GridLayoutImpl, dim, interpOrder>::deriv1D_(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
    Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const
{
    uint32 iDirX = static_cast<uint32>(Direction::X);

    uint32 iQtyOperand                = static_cast<uint32>(operand.physicalQuantity());
    auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

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

#endif // GridLayout_H
