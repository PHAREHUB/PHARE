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
 * @brief A GridLayout object represents a uniform cartesian mesh.
 * It provides methods to manipulate physical quantities discretized on such
 * a mesh. From the client code, GridLayout is the object that knows
 * the centering of all quantities in the simulation. It knows where are
 * the physical domain indexes VS the ghost indexes for a grid.
 *
 * The specific knowledge of which quantity is where is encapsulated in the private
 * GridLayout implementation GridLayoutImpl. For instance, the Yee lattice is
 * the GridLayoutImplYee of such implementations.
 *
 * GridLayout in itself gathers many methods to manipulate discretized quantities.
 * All these methods only depend on the concept of primal and dual nodes.
 * Primal nodes are nodes at integer multiple of the meshSize, dual node are positionned
 * at half multiples of the mesh size. In, say, 2D, a given quantity can for instance
 * be (primal,dual), meaning that it will be found at cell corners in the X direction
 * , but offset by dy/2 in the Y direction.
 *
 * GridLayout is templated by a specific implementation, such as GridLayoutImplYee.
 * This implementation brings the dimensionality and the interpolation order compile-time
 * constants to the GridLayout.
 *
 */
template<typename GridLayoutImpl>
class GridLayout
{
public:
    static constexpr std::size_t dimension    = GridLayoutImpl::dimension;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;



    /**
     * @brief Constructor of a GridLayout
     * @param meshSize is the size of the mesh in each direction
     * @param nbrCells is the number of physical cells of the grid
     * @param origin is the point of origin in physical units of the origin of the grid
     */
    GridLayout(std::array<double, dimension> const& meshSize,
               std::array<uint32, dimension> const& nbrCells,
               Point<double, dimension> const& origin)
        : meshSize_{meshSize}
        , origin_{origin}
        , nbrPhysicalCells_{nbrCells}
        , physicalStartIndexTable_{initPhysicalStart_()}
        , physicalEndIndexTable_{initPhysicalEnd_()}
        , ghostEndIndexTable_{initGhostEnd_()}
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




    /**
     * @brief origin return the lower point of the grid described by the GridLayout
     * in physical coordinates
     */
    Point<double, dimension> origin() const noexcept { return origin_; }




    /**
     * @brief returns the mesh size in the 'dim' dimensions
     */
    std::array<double, dimension> meshSize() const noexcept { return meshSize_; }




    double inverseMeshSize(Direction direction) const noexcept
    {
        return inverseMeshSize_[static_cast<uint32>(direction)];
    }




    std::array<double, dimension> inverseMeshSize() const noexcept { return inverseMeshSize_; }




    /**
     * @brief nbrCells returns the number of cells in the physical domain
     * described by the gridlayout
     */
    std::array<uint32, dimension> nbrCells() const { return nbrPhysicalCells_; }




    /**
     * @brief physicalStartIndex returns the index of the first node of a given
     * centering and in a given direction that is in the physical domain, i.e. not a ghost node.
     */
    uint32 physicalStartIndex(QtyCentering centering, Direction direction) const
    {
        uint32 icentering = static_cast<uint32>(centering);
        uint32 iDir       = static_cast<uint32>(direction);
        return physicalStartIndexTable_[icentering][iDir];
    }




    uint32 physicalStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                              Direction direction) const
    {
        uint32 iQty                        = static_cast<uint32>(hybridQuantity);
        uint32 iDir                        = static_cast<uint32>(direction);
        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
        uint32 iCentering                  = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);

        return physicalStartIndexTable_[iCentering][iDir];
    }




    template<typename NdArrayImpl>
    uint32 physicalStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                              Direction direction) const
    {
        return physicalStartIndex(field.physicalQuantity(), direction);
    }



    /**
     * @brief physicalEndIndex returns the index of the last node of a given
     * centering and in a given direction that is in the physical domain, i.e. not a ghost node.
     */
    uint32 physicalEndIndex(QtyCentering centering, Direction direction) const
    {
        uint32 icentering = static_cast<uint32>(centering);
        uint32 iDir       = static_cast<uint32>(direction);

        return physicalEndIndexTable_[icentering][iDir];
    }




    uint32 physicalEndIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
    {
        uint32 iQty                        = static_cast<uint32>(hybridQuantity);
        uint32 iDir                        = static_cast<uint32>(direction);
        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
        uint32 iCentering                  = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);

        return physicalEndIndexTable_[iCentering][iDir];
    }




    template<typename NdArrayImpl>
    uint32 physicalEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                            Direction direction) const
    {
        return physicalEndIndex(field.physicalQuantity(), direction);
    }



    /**
     * @brief ghostStartIndex retuns the index of the first ghost node of a given centering
     * in a given direction. This is always zero by convention. This function exists only
     * for readibility reasons, not to have literal '0' in the code.
     */
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




    /**
     * @brief ghostEndIndex returns the index of the last ghost node of a given centering
     * and in a given direction.
     */
    uint32 ghostEndIndex(QtyCentering centering, Direction direction) const
    {
        uint32 iCentering = static_cast<uint32>(centering);
        uint32 iDir       = static_cast<uint32>(direction);

        return ghostEndIndexTable_[iCentering][iDir];
    }




    uint32 ghostEndIndex(HybridQuantity::Scalar const& hybridQuantity, Direction direction) const
    {
        uint32 iQty                        = static_cast<uint32>(hybridQuantity);
        uint32 iDir                        = static_cast<uint32>(direction);
        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
        uint32 iCentering                  = static_cast<uint32>(hybridQtyCentering[iQty][iDir]);
        return ghostEndIndexTable_[iCentering][iDir];
    }




    template<typename NdArrayImpl>
    uint32 ghostEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                         Direction direction) const
    {
        return ghostEndIndex(field.physicalQuantity(), direction);
    }




    /**
     * @brief fieldNodeCoordinates returns the coordinate of a multidimensional index
     * associated with a given Field, in physical coordinates.
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


        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

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



    /**
     * @brief cellCenteredCoordinates returns the coordinates in physical units
     * of a multidimensional index that is cell-centered.
     */
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
     * @brief the number of ghost nodes on each side of the mesh for a given centering
     */
    uint32 static nbrGhosts(QtyCentering centering)
    {
        uint32 nbrGhosts = nbrPrimalGhosts_();

        if (centering == QtyCentering::dual)
        {
            nbrGhosts = nbrDualGhosts_();
        }
        return nbrGhosts;
    }



    /**
     * @brief changeCentering changes primal into dual and vice versa.
     */
    QtyCentering changeCentering(QtyCentering centering) const
    {
        QtyCentering newCentering = QtyCentering::primal;

        if (centering == QtyCentering::primal)
        {
            newCentering = QtyCentering::dual;
        }

        return newCentering;
    }



    /**
     * @brief nextIndex returns the index of the next node of a given centering
     * from an index of the opposite centering.
     * For instance, if 'centering' is Qtycentering::dual, then 'indexCenter'
     * is interpreted as a primal index and the function returns the index
     * of the dual node just on the right of that primal node. This function is useful
     * for calculating derivatives.
     * The next index is not just indexCenter+1 because this depends on the number
     * of ghost nodes for dual and primal nodes.
     */
    auto static nextIndex(QtyCentering centering, uint32 indexCenter)
    {
        return indexCenter + nextIndexTable_[centering2int(centering)];
    }



    /**
     * @brief prevIndex does the same thing as nextIndex but returns the index
     * of the node of a given centering just to the left of indexCenter.
     */
    auto static prevIndex(QtyCentering centering, uint32 indexCenter)
    {
        return indexCenter + prevIndexTable_[centering2int(centering)];
    }



    /** @brief returns the local 1st order derivative of the Field operand
     * at a multidimensional index and in a given direction.
     * The function can perform 1D, 2D and 3D 1st order derivatives, depending
     * on the dimensionality of the GridLayout.
     */
    template<typename Field, typename DirectionTag>
    auto deriv(Field const& operand, MeshIndex<Field::dimension> index, DirectionTag)
    {
        auto fieldCentering   = centering(operand.physicalQuantity());
        constexpr uint32 dirX = 0, dirY = 1, dirZ = 2;

        if constexpr (Field::dimension == 1)
        {
            auto ip  = nextIndex(fieldCentering[dirX], index.i);
            auto im  = prevIndex(fieldCentering[dirX], index.i);
            auto der = inverseMeshSize_[dirX]
                       * (operand(nextIndex(fieldCentering[dirX], index.i))
                          - operand(prevIndex(fieldCentering[dirX], index.i)));
            return der;
        }

        else if constexpr (Field::dimension == 2)
        {
            if constexpr (DirectionTag::direction == Direction::X)
            {
                auto next = operand(nextIndex(fieldCentering[dirX], index.i), index.j);
                auto prev = operand(prevIndex(fieldCentering[dirX], index.i), index.j);
                return inverseMeshSize_[dirX] * (next - prev);
            }

            if constexpr (DirectionTag::direction == Direction::Y)
            {
                auto next = operand(index.i, nextIndex(fieldCentering[dirY], index.j));
                auto prev = operand(index.i, prevIndex(fieldCentering[dirY], index.j));
                return inverseMeshSize_[dirY] * (next - prev);
            }
        }
        else if constexpr (Field::dimension == 3)
        {
            if constexpr (DirectionTag::direction == Direction::X)
            {
                auto next = operand(nextIndex(fieldCentering[dirX], index.i), index.j, index.k);
                auto prev = operand(prevIndex(fieldCentering[dirX], index.i), index.j, index.k);
                return inverseMeshSize_[dirX] * (next - prev);
            }

            if constexpr (DirectionTag::direction == Direction::Y)
            {
                auto next = operand(index.i, nextIndex(fieldCentering[dirY], index.j), index.k);
                auto prev = operand(index.i, prevIndex(fieldCentering[dirY], index.j), index.k);
                return inverseMeshSize_[dirY] * (next - prev);
            }
            if constexpr (DirectionTag::direction == Direction::Z)
            {
                auto next = operand(index.i, index.j, nextIndex(fieldCentering[dirZ], index.k));
                auto prev = operand(index.i, index.j, prevIndex(fieldCentering[dirZ], index.k));
                return inverseMeshSize_[dirZ] * (next - prev);
            }
        }
    }



    // ----------------------------------------------------------------------
    //                      LAYOUT SPECIFIC METHODS
    //
    // The methods below return results that need information specific to the
    // layout that is used. They thus all refer to the GridLayoutImpl.
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

        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        std::array<QtyCentering, dimension> qtyCentering;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
        {
            qtyCentering[iDir] = hybridQtyCentering[iQty][iDir];
        }

        return nodeNbrFromCentering_(qtyCentering);
    }



    /**
     * @brief allocSizeDerived returns the shape of the array to be allocated to store
     * the derivative of a given quantity in a given direction.
     */
    std::array<uint32, dimension> allocSizeDerived(HybridQuantity::Scalar qty, Direction dir) const
    {
        uint32 iDerivedDir = static_cast<uint32>(dir);
        uint32 iQty        = static_cast<uint32>(qty);

        // get the centering of the derivative of 'qty' in the direction of derivation
        QtyCentering newCentering = derivedCentering(qty, dir);

        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

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




    /** @brief return the centering of a given Field along a given direction
     */
    template<typename NdArrayImpl>
    QtyCentering fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                Direction dir) const
    {
        uint32 iDir                        = static_cast<uint32>(dir);
        uint32 iQty                        = static_cast<uint32>(field.physicalQuantity());
        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        return hybridQtyCentering[iQty][iDir];
    }



    /**
     * @brief nbrPhysicalNodes returns the number of nodes in each direction, that are node ghost
     * nodes
     */
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
     * @brief derivedCentering this function returns the
     * centering (primal or dual) of a quantity after a first order derivation. dual becomes primal
     * and primal becomes dual. hybridQuantityCentering is used to know if the
     * HybridQuantity::Quantity 'qty' is primal or dual in the Direction 'dir'
     */
    QtyCentering derivedCentering(HybridQuantity::Scalar qty, Direction dir) const
    {
        uint32 iField = static_cast<uint32>(qty);
        uint32 idir   = static_cast<uint32>(dir);


        constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

        QtyCentering newCentering = changeCentering(hybridQtyCentering[iField][idir]);

        return newCentering;
    }



    /**
     * @brief momentsToEx return the indexes and associated coef to compute the linear
     * interpolation necessary to project moments onto Ex.
     */
    auto static constexpr momentsToEx() { return GridLayoutImpl::momentsToEx(); }



    /**
     * @brief momentsToEy return the indexes and associated coef to compute the linear
     * interpolation necessary to project moments onto Ey.
     */
    auto static constexpr momentsToEy() { return GridLayoutImpl::momentsToEy(); }



    /**
     * @brief momentsToEz return the indexes and associated coef to compute the linear
     * interpolation necessary to project moments onto Ez.
     */
    auto static constexpr momentsToEz() { return GridLayoutImpl::momentsToEz(); }




    /**
     * @brief ExToMoments return the indexes and associated coef to compute the linear
     * interpolation necessary to project Ex onto moments.
     */
    auto static constexpr ExToMoments() { return GridLayoutImpl::ExToMoments(); }




    /**
     * @brief EyToMoments return the indexes and associated coef to compute the linear
     * interpolation necessary to project Ey onto moments.
     */
    auto static constexpr EyToMoments() { return GridLayoutImpl::EyToMoments(); }




    /**
     * @brief EzToMoments return the indexes and associated coef to compute the linear
     * interpolation necessary to project Ez onto moments.
     */
    auto static constexpr EzToMoments() { return GridLayoutImpl::EzToMoments(); }




    /**
     * @brief ByToEx return the indexes and associated coef to compute the linear
     * interpolation necessary to project By onto Ex.
     */
    auto static constexpr ByToEx() { return GridLayoutImpl::ByToEx(); }




    /**
     * @brief BzToEx return the indexes and associated coef to compute the linear
     * interpolation necessary to project Bz onto Ex.
     */
    auto static constexpr BzToEx() { return GridLayoutImpl::BzToEx(); }




    /**
     * @brief BxToEy return the indexes and associated coef to compute the linear
     * interpolation necessary to project Bx onto Ey.
     */
    auto static constexpr BxToEy() { return GridLayoutImpl::BxToEy(); }




    /**
     * @brief BzToEy return the indexes and associated coef to compute the linear
     * interpolation necessary to project Bz onto Ey.
     */
    auto static constexpr BzToEy() { return GridLayoutImpl::BzToEy(); }




    /**
     * @brief BxToEz return the indexes and associated coef to compute the linear
     * interpolation necessary to project Bx onto Ez.
     */
    auto static constexpr BxToEz() { return GridLayoutImpl::BxToEz(); }




    /**
     * @brief ByToEz return the indexes and associated coef to compute the linear
     * interpolation necessary to project By onto Ez.
     */
    auto static constexpr ByToEz() { return GridLayoutImpl::ByToEz(); }


private:
    /**
     * @brief nextPrimal_ returns the index shift needed to go to the next primal
     * node from a dual node. This depends on whether the dual have more ghost nodes
     * than the primal.
     * returning 0 means that the next primal has the same index as the current dual.
     * returning 1 means that the next primal index is the current dual + 1
     */
    constexpr static auto nextPrimal_()
    {
        if constexpr (nbrDualGhosts_() > nbrPrimalGhosts_())
        {
            return 0;
        }
        else if constexpr (nbrDualGhosts_() == nbrPrimalGhosts_())
        {
            return 1;
        }
    }



    /**
     * @brief prevPrimal_ does the same as nextPrimal_ but for the previous primal
     */
    constexpr static auto prevPrimal_()
    {
        if constexpr (nbrDualGhosts_() > nbrPrimalGhosts_())
        {
            return -1;
        }
        else if constexpr (nbrDualGhosts_() == nbrPrimalGhosts_())
        {
            return 0;
        }
    }



    /**
     * @brief nextDual_ is identical to nextPrimal for dual nodes
     */
    constexpr static auto nextDual_()
    {
        if constexpr (nbrDualGhosts_() > nbrPrimalGhosts_())
        {
            return 1;
        }
        else if constexpr (nbrDualGhosts_() == nbrPrimalGhosts_())
        {
            return 0;
        }
    }



    /**
     * @brief prevDual_ is identical to prevPrimal_ for dual nodes.
     */
    constexpr static auto prevDual_()
    {
        if constexpr (nbrDualGhosts_() > nbrPrimalGhosts_())
        {
            return 0;
        }
        else if constexpr (nbrDualGhosts_() == nbrPrimalGhosts_())
        {
            return -1;
        }
    }


    /**
     * @brief nbrDualGhosts_ returns the number of ghost nodes on each side for dual quantities.
     * The formula is based only on the interpolation order, whch means only particle-mesh
     * interactions constrain the number of dual ghost nodes.
     */
    auto constexpr static nbrDualGhosts_() { return static_cast<uint32>((interp_order + 1) / 2.); }


    /**
     * @brief nbrPrimalGhosts_ returns the number of primal ghost nodes.
     * Contrary to dual ghost nodes, the formula to get the number of primal ghost nodes depend on
     * the interpolation order. If based only on the particle-mesh interaction, order1 would
     * not need primal ghost nodes. But there is a minimum of 1 ghost node for primal that is linked
     * to the possibility of calculating second order derivatives of primal quantities (e.g.
     * laplacian of J for a yee lattice). Dual ghosts don't have this issue since they always have
     * at least 1 ghost.
     */
    auto constexpr static nbrPrimalGhosts_()
    {
        if constexpr (interp_order == 1)
        {
            return nbrDualGhosts_();
        }
        else if constexpr (interp_order == 2 || interp_order == 3)
        {
            return static_cast<uint32>(interp_order / 2.);
        }
    }




    uint32 static constexpr dualOffset_() noexcept { return 1; }



    /**
     * @brief physicalNodeNbrFromCentering_ returns the number of physical nodes for all directions
     * depending on the multi-dimensional centering.
     */
    std::array<uint32, dimension>
    physicalNodeNbrFromCentering_(std::array<QtyCentering, dimension> const& qtyCenterings) const
    {
        std::array<uint32, dimension> nodeNbr;

        for (std::size_t iDir = 0; iDir < dimension; ++iDir)
        {
            nodeNbr[iDir] = nbrPhysicalCells_[iDir] + 1
                            - ((qtyCenterings[iDir] == QtyCentering::dual) ? dualOffset_() : 0);
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



    auto initPhysicalStart_()
    {
        std::array<std::array<uint32, dimension>, 2> physicalStartIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        physicalStartIndexTable[iprimal][data.idirX] = nbrPrimalGhosts_();
        physicalStartIndexTable[idual][data.idirX]   = nbrDualGhosts_();

        if constexpr (dimension > 1)
        {
            physicalStartIndexTable[iprimal][data.idirY] = nbrPrimalGhosts_();
            physicalStartIndexTable[idual][data.idirY]   = nbrDualGhosts_();

            if constexpr (dimension > 2)
            {
                physicalStartIndexTable[iprimal][data.idirZ] = nbrPrimalGhosts_();
                physicalStartIndexTable[idual][data.idirZ]   = nbrDualGhosts_();
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
    auto initPhysicalEnd_()
    {
        std::array<std::array<uint32, dimension>, 2> physicalEndIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        physicalEndIndexTable[iprimal][data.idirX]
            = physicalStartIndexTable_[iprimal][data.idirX] + nbrPhysicalCells_[data.idirX];


        physicalEndIndexTable[idual][data.idirX] = physicalStartIndexTable_[idual][data.idirX]
                                                   + nbrPhysicalCells_[data.idirX] - dualOffset_();

        if constexpr (dimension > 1)
        {
            physicalEndIndexTable[iprimal][data.idirY]
                = physicalStartIndexTable_[iprimal][data.idirY] + nbrPhysicalCells_[data.idirY];

            physicalEndIndexTable[idual][data.idirY] = physicalStartIndexTable_[idual][data.idirY]
                                                       + nbrPhysicalCells_[data.idirY]
                                                       - dualOffset_();

            if constexpr (dimension > 2)
            {
                physicalEndIndexTable[iprimal][data.idirZ]
                    = physicalStartIndexTable_[iprimal][data.idirZ] + nbrPhysicalCells_[data.idirZ];

                physicalEndIndexTable[idual][data.idirZ]
                    = physicalStartIndexTable_[idual][data.idirZ] + nbrPhysicalCells_[data.idirZ]
                      - dualOffset_();
            }
        }
        return physicalEndIndexTable;
    }




    /**
     * @brief GridLayout<GridLayoutImpl::dim>::initGhostEnd calculate and stores the index
     * of the last primal and dual nodes in each direction. The formula simply
     * consists in starting at physicalEndIndex() and to add the number of ghost nodes.
     */
    auto initGhostEnd_()
    {
        std::array<std::array<uint32, dimension>, 2> ghostEndIndexTable;

        uint32 iprimal = static_cast<uint32>(data.primal);
        uint32 idual   = static_cast<uint32>(data.dual);

        ghostEndIndexTable[iprimal][data.idirX]
            = physicalEndIndexTable_[iprimal][data.idirX] + nbrPrimalGhosts_();

        ghostEndIndexTable[idual][data.idirX]
            = physicalEndIndexTable_[idual][data.idirX] + nbrDualGhosts_();

        if constexpr (dimension > 1)
        {
            ghostEndIndexTable[iprimal][data.idirY]
                = physicalEndIndexTable_[iprimal][data.idirY] + nbrPrimalGhosts_();

            ghostEndIndexTable[idual][data.idirY]
                = physicalEndIndexTable_[idual][data.idirY] + nbrDualGhosts_();

            if constexpr (dimension > 2)
            {
                ghostEndIndexTable[iprimal][data.idirZ]
                    = physicalEndIndexTable_[iprimal][data.idirZ] + nbrPrimalGhosts_();

                ghostEndIndexTable[idual][data.idirZ]
                    = physicalEndIndexTable_[idual][data.idirZ] + nbrDualGhosts_();
            }
        }
        return ghostEndIndexTable;
    }




    std::array<double, dimension> meshSize_;
    Point<double, dimension> origin_;
    std::array<uint32, dimension> nbrPhysicalCells_;
    std::array<double, dimension> inverseMeshSize_;
    static constexpr gridDataT data{};

    // stores key indices in each direction (3) for primal and dual nodes (2)
    std::array<std::array<uint32, dimension>, 2> physicalStartIndexTable_;
    std::array<std::array<uint32, dimension>, 2> physicalEndIndexTable_;
    std::array<std::array<uint32, dimension>, 2> ghostEndIndexTable_;


    // this constexpr initialization only works if primal==0 and dual==1
    // this is defined in gridlayoutdefs.h don't change it because these
    // arrays will be accessed with [primal] and [dual] indexes.
    constexpr static std::array<int, 2> nextIndexTable_{{nextPrimal_(), nextDual_()}};
    constexpr static std::array<int, 2> prevIndexTable_{{prevPrimal_(), prevDual_()}};
};




} // namespace PHARE

#endif // GridLayout_H
