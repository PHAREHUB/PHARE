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
#include "utilities/index/index.h"
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
template<typename GridLayoutImpl>
class GridLayout
{
    // ------------------------------------------------------------------------
    //                             PROTECTED
    //
    // this code will be shared by all concrete of GridLayoutImpl*
    // ------------------------------------------------------------------------
private:
    static constexpr std::size_t dimension    = GridLayoutImpl::dimension;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;

public:
    /* uint32 static constexpr nbdims_{GridLayoutImpl::getDimensions()}; */
    // uint32 constexpr static nbdims_{dim};
    std::array<double, dimension> meshSize_;
    Point<double, dimension> origin_;
    std::array<uint32, dimension> nbrPhysicalCells_;

    static constexpr gridDataT data{};



    // stores key indices in each direction (3) for primal and dual nodes (2)
    std::array<std::array<uint32, dimension>, 2> physicalStartIndexTable_;
    std::array<std::array<uint32, dimension>, 2> physicalEndIndexTable_;
    std::array<std::array<uint32, dimension>, 2> ghostEndIndexTable_;

    std::array<double, dimension> inverseMeshSize_;

    GridLayout(std::array<double, dimension> const& meshSize,
               std::array<uint32, dimension> const& nbrCells,
               Point<double, dimension> const& origin)
        : meshSize_{meshSize}
        , origin_{origin}
        , nbrPhysicalCells_{nbrCells}
        , physicalStartIndexTable_{initPhysicalStart()}
        , physicalEndIndexTable_{initPhysicalEnd()}
        , ghostEndIndexTable_{initGhostEnd()}
    {
        inverseMeshSize_[0] = 1. / meshSize_[0];
        if constexpr (dimension > 1)
        {
            inverseMeshSize_[1] = 1. / meshSize_[1];
            if constexpr (dimension > 2)
            {
                inverseMeshSize_[2] = 1. / meshSize_[2];
            }
        }
    }

    auto initPhysicalStart()
    {
        std::array<std::array<uint32, dimension>, 2> physicalStartIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        physicalStartIndexTable[iprimal][data.idirX] = nbrPrimalGhosts();
        physicalStartIndexTable[idual][data.idirX]   = nbrDualGhosts();

        if constexpr (dimension > 1)
        {
            physicalStartIndexTable[iprimal][data.idirY] = nbrPrimalGhosts();
            physicalStartIndexTable[idual][data.idirY]   = nbrDualGhosts();

            if constexpr (dimension > 2)
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
        std::array<std::array<uint32, dimension>, 2> physicalEndIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        physicalEndIndexTable[iprimal][data.idirX]
            = physicalStartIndexTable_[iprimal][data.idirX] + nbrPhysicalCells_[data.idirX];


        physicalEndIndexTable[idual][data.idirX] = physicalStartIndexTable_[idual][data.idirX]
                                                   + nbrPhysicalCells_[data.idirX] - dualOffset();

        if constexpr (dimension > 1)
        {
            physicalEndIndexTable[iprimal][data.idirY]
                = physicalStartIndexTable_[iprimal][data.idirY] + nbrPhysicalCells_[data.idirY];

            physicalEndIndexTable[idual][data.idirY] = physicalStartIndexTable_[idual][data.idirY]
                                                       + nbrPhysicalCells_[data.idirY]
                                                       - dualOffset();

            if constexpr (dimension > 2)
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
        std::array<std::array<uint32, dimension>, 2> ghostEndIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        ghostEndIndexTable[iprimal][data.idirX]
            = physicalEndIndexTable_[iprimal][data.idirX] + nbrPrimalGhosts();

        ghostEndIndexTable[idual][data.idirX]
            = physicalEndIndexTable_[idual][data.idirX] + nbrDualGhosts();

        if constexpr (dimension > 1)
        {
            ghostEndIndexTable[iprimal][data.idirY]
                = physicalEndIndexTable_[iprimal][data.idirY] + nbrPrimalGhosts();

            ghostEndIndexTable[idual][data.idirY]
                = physicalEndIndexTable_[idual][data.idirY] + nbrDualGhosts();

            if constexpr (dimension > 2)
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
    Point<double, dimension> origin() const { return origin_; }


    /**
     * @brief returns the mesh size in the 'dim' dimensions
     */
    std::array<double, dimension> meshSize() const { return meshSize_; }


    double inverseMeshSize(Direction direction) const noexcept
    {
        return inverseMeshSize_[static_cast<uint32>(direction)];
    }


    std::array<double, dimension> inverseMeshSize() const noexcept { return inverseMeshSize_; }


    std::array<uint32, dimension> nbrCells() const { return nbrPhysicalCells_; }



    uint32 nbDimensions() const { return dimension; }



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
    template<std::size_t order                         = interp_order,
             std::enable_if_t<order == 1, dummy::type> = dummy::value>
    uint32 constexpr nbrDualGhosts() const
    {
        /* for first order Interpolation, there is no primal ghost node neeeded
           for particle/mesh interactions. However one ghost node is required
           for calculating Laplacians so we add one.
        */
        return static_cast<uint32>((interp_order + 1) / 2.);
    }


    template<std::size_t order                         = interp_order,
             std::enable_if_t<order == 1, dummy::type> = dummy::value>
    uint32 constexpr nbrPrimalGhosts() const
    {
        /* for first order Interpolation, there is no primal ghost node neeeded
           for particle/mesh interactions. However one ghost node is required
           for calculating Laplacians so we add one.
        */
        return static_cast<uint32>((interp_order + 1) / 2.);
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
    template<std::size_t order                          = interp_order,
             std::enable_if_t<(order > 1), dummy::type> = dummy::value>
    uint32 constexpr nbrDualGhosts() const
    {
        /* for interpolation order larger than 1, there is at least 1 primal ghost
           node so Laplacians can be calculated OK */
        return static_cast<uint32>((interp_order + 1) / 2.);
    }

    template<std::size_t order                          = interp_order,
             std::enable_if_t<(order > 1), dummy::type> = dummy::value>
    uint32 constexpr nbrPrimalGhosts() const
    {
        /* for interpolation order larger than 1, there is at least 1 primal ghost
           node so Laplacians can be calculated OK */
        return static_cast<uint32>(interp_order / 2.);
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
    Point<double, dimension>
    fieldNodeCoordinates(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                         const Point<double, dimension>& origin, Indexes... index) const
    {
        static_assert(sizeof...(Indexes) == dimension,
                      "Error dimension does not match number of arguments");


        uint32 iQuantity       = static_cast<uint32>(field.physicalQuantity());
        constexpr uint32 iDual = static_cast<uint32>(QtyCentering::dual);


        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        Point<int32, dimension> coord{static_cast<int32>(index)...};

        Point<double, dimension> position;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
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
    Point<double, dimension> cellCenteredCoordinates(Indexes... index) const
    {
        static_assert(sizeof...(Indexes) == dimension,
                      "Error dimension does not match number of arguments");

        uint32 constexpr iPrimal = static_cast<uint32>(QtyCentering::primal);

        constexpr double halfCell = 0.5;
        // A shift of +dx/2, +dy/2, +dz/2 is necessary to get the
        // cell center physical coordinates,
        // because this point is located on the dual mesh

        Point<uint32, dimension> coord(index...);

        Point<double, dimension> physicalPosition;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
        {
            auto iStart = physicalStartIndexTable_[iPrimal][iDir];

            physicalPosition[iDir]
                = (static_cast<double>(coord[iDir] - iStart) + halfCell) * meshSize_[iDir]
                  + origin_[iDir];
        }

        return physicalPosition;
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
    constexpr static std::array<QtyCentering, dimension>
    centering(HybridQuantity::Scalar const& hybridQuantity)
    {
        return GridLayoutImpl::centering(hybridQuantity);
    }




    /**
     * @brief GridLayout<GridLayoutImpl::dim>::allocSize_
     * @return An std::array<uint32, dim> object, containing the size to which allocate arrays
     * of an HybridQuantity::Quantity 'qty' in every directions.
     */
    std::array<uint32, dimension> allocSize(HybridQuantity::Scalar qty) const
    {
        uint32 iQty = static_cast<uint32>(qty);


        // TODO: hybridQtyCentering should be defined per dimension so that we could simply do
        // auto sizeArray = nodeNbrFromCentering_(hybridQtyCentering[iQty]);

        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        std::array<QtyCentering, dimension> qtyCentering;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
        {
            qtyCentering[iDir] = hybridQtyCentering[iQty][iDir];
        }

        return nodeNbrFromCentering_(qtyCentering);
    }



    std::array<uint32, dimension> allocSizeDerived(HybridQuantity::Scalar qty, Direction dir) const
    {
        uint32 iDerivedDir = static_cast<uint32>(dir);
        uint32 iQty        = static_cast<uint32>(qty);

        // get the centering of the derivative of 'qty' in the direction of derivation
        QtyCentering newCentering = derivedCentering(qty, dir);

        auto constexpr hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        std::array<QtyCentering, dimension> qtyCenterings;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
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




    std::array<uint32, dimension> nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const
    {
        std::array<QtyCentering, dimension> centerings;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
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


    constexpr auto nextPrimal(int dualIndex)
    {
        if constexpr (nbrDualGhosts() > nbrPrimalGhosts())
        {
            return dualIndex + 1;
        }
        else if constexpr (nbrDualGhosts() == nbrPrimalGhosts())
        {
            return dualIndex;
        }
    }

    constexpr auto prevPrimal(int dualIndex)
    {
        if constexpr (nbrDualGhosts() > nbrPrimalGhosts())
        {
            return dualIndex;
        }
        else if constexpr (nbrDualGhosts() == nbrPrimalGhosts())
        {
            return dualIndex - 1;
        }
    }



    constexpr auto nextDual(int primalIndex)
    {
        if constexpr (nbrDualGhosts() > nbrPrimalGhosts())
        {
            return primalIndex;
        }
        else if constexpr (nbrDualGhosts() == nbrPrimalGhosts())
        {
            return primalIndex + 1;
        }
    }


    constexpr auto prevDual(int primalIndex)
    {
        if constexpr (nbrDualGhosts() > nbrPrimalGhosts())
        {
            return primalIndex;
        }
        else if constexpr (nbrDualGhosts() == nbrPrimalGhosts())
        {
            return primalIndex - 1;
        }
    }


    template<typename Field, typename DirectionTag>
    auto deriv(Field const& operand, MeshIndex<Field::dimension> index, DirectionTag)
    {
        constexpr auto fieldCentering = centering(operand.physicalQuantity());

        if constexpr (Field::dimension == 1)
        {
            if constexpr (fieldCentering[0] == QtyCentering::dual)
            {
                return inverseMeshSize_[0] * operand(nextPrimal(index.i))
                       - operand(prevPrimal(index.i));
            }

            else if constexpr (fieldCentering[0] == QtyCentering::primal)
            {
                return inverseMeshSize_[0] * operand(nextDual(index.i))
                       - operand(prevDual(index.i));
            }
        }
        if constexpr (Field::dimension == 2)
        {
            if constexpr (DirectionTag::direction == Direction::X)
            {
                if constexpr (fieldCentering[0] == QtyCentering::dual)
                {
                    return inverseMeshSize_[0] * operand(index.i, nextPrimal(index.j))
                           - operand(prevPrimal(index.j));
                }

                else if constexpr (fieldCentering[0] == QtyCentering::primal)
                {
                    return inverseMeshSize_[0] * operand(index.i, nextDual(index.j))
                           - operand(index.i, prevDual(index.j));
                }
            }
            else if constexpr (DirectionTag::direction == Direction::Y)
            {
                if constexpr (fieldCentering[1] == QtyCentering::dual)
                {
                    return inverseMeshSize_[1] * operand(index.i, nextPrimal(index.j))
                           - operand(index.i, prevPrimal(index.j));
                }

                else if constexpr (fieldCentering[1] == QtyCentering::primal)
                {
                    return inverseMeshSize_[1] * operand(index.i, nextDual(index.j))
                           - operand(prevDual(index.j));
                }
            }
        }


        if constexpr (Field::dimension == 3)
        {
            if constexpr (DirectionTag::direction == Direction::X)
            {
                if constexpr (fieldCentering[0] == QtyCentering::dual)
                {
                    return inverseMeshSize_[0] * operand(nextPrimal(index.i), index.j, index.k)
                           - operand(prevPrimal(index.i), index.j, index.k);
                }

                else if constexpr (fieldCentering[0] == QtyCentering::primal)
                {
                    return inverseMeshSize_[0] * operand(nextDual(index.i), index.j, index.k)
                           - operand(prevDual(index.i), index.j, index.k);
                }
            }
            else if constexpr (DirectionTag::direction == Direction::Y)
            {
                if constexpr (fieldCentering[1] == QtyCentering::dual)
                {
                    return inverseMeshSize_[1] * operand(index.i, nextPrimal(index.j), index.k)
                           - operand(index.i, prevPrimal(index.j), index.k);
                }

                else if constexpr (fieldCentering[1] == QtyCentering::primal)
                {
                    return inverseMeshSize_[1] * operand(index.i, nextDual(index.j), index.k)
                           - operand(index.i, prevDual(index.j), index.k);
                }
            }

            else if constexpr (DirectionTag::direction == Direction::Z)
            {
                if constexpr (fieldCentering[2] == QtyCentering::dual)
                {
                    return inverseMeshSize_[2] * operand(index.i, index.j, nextPrimal(index.k))
                           - operand(index.i, index.j, prevPrimal(index.k));
                }

                else if constexpr (fieldCentering[2] == QtyCentering::primal)
                {
                    return inverseMeshSize_[2] * operand(index.i, index.j, nextDual(index.k))
                           - operand(index.i, index.j, prevDual(index.i));
                }
            } // 3D directionZ
        }     // 3D
    }



    auto static constexpr momentsToEx() { return GridLayoutImpl::momentsToEx(); }
    auto static constexpr momentsToEy() { return GridLayoutImpl::momentsToEy(); }
    auto static constexpr momentsToEz() { return GridLayoutImpl::momentsToEz(); }
    auto static constexpr ExToMoments() { return GridLayoutImpl::ExToMoments(); }
    auto static constexpr EyToMoments() { return GridLayoutImpl::EyToMoments(); }
    auto static constexpr EzToMoments() { return GridLayoutImpl::EzToMoments(); }
    auto static constexpr ByToEx() { return GridLayoutImpl::ByToEx(); }
    auto static constexpr BzToEx() { return GridLayoutImpl::BzToEx(); }
    auto static constexpr BxToEy() { return GridLayoutImpl::BxToEy(); }
    auto static constexpr BzToEy() { return GridLayoutImpl::BzToEy(); }
    auto static constexpr BxToEz() { return GridLayoutImpl::BxToEz(); }
    auto static constexpr ByToEz() { return GridLayoutImpl::ByToEz(); }


private:
    std::array<uint32, dimension>
    physicalNodeNbrFromCentering_(std::array<QtyCentering, dimension> const& qtyCenterings) const
    {
        std::array<uint32, dimension> nodeNbr;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
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
    std::array<uint32, dimension>
    nodeNbrFromCentering_(std::array<QtyCentering, dimension> const& qtyCenterings) const
    {
        std::array<uint32, dimension> nbrNodes = physicalNodeNbrFromCentering_(qtyCenterings);

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
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
#if 0
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
#endif



} // namespace PHARE

#endif // GridLayout_H
