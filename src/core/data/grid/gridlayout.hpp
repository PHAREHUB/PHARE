#ifndef PHARE_CORE_GRID_GridLayout_HPP
#define PHARE_CORE_GRID_GridLayout_HPP


#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "gridlayoutdefs.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/point/point.hpp"
#include "core/def.hpp"

#include <array>
#include <cmath>
#include <tuple>
#include <cstddef>
#include <functional>
#include <type_traits>

namespace PHARE
{
namespace core
{
    template<typename T, typename Attempt = void>
    struct has_physicalQuantity : std::false_type
    {
    };

    template<typename T>
    struct has_physicalQuantity<
        T, core::tryToInstanciate<decltype(std::declval<T>().physicalQuantity())>> : std::true_type
    {
    };
    template<typename T>
    constexpr bool has_physicalQuantity_v = has_physicalQuantity<T>::value;


    NO_DISCARD constexpr int centering2int(QtyCentering c)
    {
        return static_cast<int>(c);
    }


    template<std::size_t interpOrder>
    NO_DISCARD std::uint32_t constexpr ghostWidthForParticles()
    {
        return (interpOrder % 2 == 0 ? interpOrder / 2 + 1 : (interpOrder + 1) / 2);
    }



    template<typename T, std::size_t s>
    NO_DISCARD auto boxFromNbrCells(std::array<T, s> nbrCells)
    {
        Point<int, s> lower;
        Point<int, s> upper;

        for (auto i = 0u; i < nbrCells.size(); ++i)
        {
            lower[i] = 0;
            upper[i] = static_cast<int>(nbrCells[i]) - 1;
        };

        return Box{lower, upper};
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
    protected:
        GridLayout() = default;

    public:
        static constexpr std::size_t dimension    = GridLayoutImpl::dimension;
        static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
        using This                                = GridLayout<GridLayoutImpl>;
        using implT                               = GridLayoutImpl;


        /**
         * @brief Constructor of a GridLayout
         * @param meshSize is the size of the mesh in each direction
         * @param nbrCells is the number of physical cells of the grid
         * @param origin is the point of origin in physical units of the origin of the grid
         */
        GridLayout(std::array<double, dimension> const& meshSize,
                   std::array<std::uint32_t, dimension> const& nbrCells,
                   Point<double, dimension> const& origin,
                   Box<int, dimension> AMRBox = Box<int, dimension>{})
            : meshSize_{meshSize}
            , origin_{origin}
            , nbrPhysicalCells_{nbrCells}
            , physicalStartIndexTable_{initPhysicalStart_()}
            , physicalEndIndexTable_{initPhysicalEnd_()}
            , ghostEndIndexTable_{initGhostEnd_()}
            , AMRBox_{AMRBox}
        {
            if (AMRBox_.isEmpty())
            {
                AMRBox_ = boxFromNbrCells(nbrCells);
            }
            else
            {
                if (AMRBox.size() != boxFromNbrCells(nbrCells).size())
                {
                    throw std::runtime_error("Error - invalid AMR box, incorrect number of cells");
                }
            }


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


        GridLayout(GridLayout const& that) = default;
        GridLayout(GridLayout&& source)    = default;


        /**
         * @brief origin return the lower point of the grid described by the GridLayout
         * in physical coordinates
         */
        NO_DISCARD Point<double, dimension> origin() const noexcept { return origin_; }



        /**
         * @brief returns the mesh size in the 'dim' dimensions
         */
        NO_DISCARD std::array<double, dimension> const& meshSize() const noexcept
        {
            return meshSize_;
        }



        NO_DISCARD double inverseMeshSize(Direction direction) const noexcept
        {
            return inverseMeshSize_[static_cast<std::uint32_t>(direction)];
        }



        NO_DISCARD std::array<double, dimension> inverseMeshSize() const noexcept
        {
            return inverseMeshSize_;
        }



        /**
         * @brief nbrCells returns the number of cells in the physical domain
         * described by the gridlayout
         */
        NO_DISCARD auto& nbrCells() const { return nbrPhysicalCells_; }


        NO_DISCARD auto const& AMRBox() const { return AMRBox_; }



        NO_DISCARD static std::size_t constexpr nbrParticleGhosts()
        {
            return ghostWidthForParticles<interp_order>();
        }

        template<typename Centering, typename Direction>
        NO_DISCARD auto ghostStartToEnd(Centering const& centering, Direction const direction) const
        {
            return std::make_tuple(ghostStartIndex(centering, direction),
                                   ghostEndIndex(centering, direction));
        }

        template<typename Field, std::enable_if_t<has_physicalQuantity_v<Field>, bool> = 0>
        NO_DISCARD auto ghostStartToEnd(Field const& field, Direction direction) const
        {
            return ghostStartToEnd(field.physicalQuantity(), direction);
        }

        template<typename Centering, typename Direction>
        NO_DISCARD auto physicalStartToEnd(Centering const& centering,
                                           Direction const direction) const
        {
            return std::make_tuple(physicalStartIndex(centering, direction),
                                   physicalEndIndex(centering, direction));
        }


        template<typename Field, std::enable_if_t<has_physicalQuantity_v<Field>, bool> = 0>
        NO_DISCARD auto physicalStartToEnd(Field const& field, Direction direction) const
        {
            return physicalStartToEnd(field.physicalQuantity(), direction);
        }



        template<typename Centering>
        NO_DISCARD auto physicalStartToEndIndices(Centering const& centering,
                                                  bool const includeEnd = false) const
        {
            return StartToEndIndices_(
                centering,
                [&](auto const& centering_, auto const direction) {
                    return this->physicalStartToEnd(centering_, direction);
                },
                includeEnd);
        }

        template<typename Centering>
        NO_DISCARD auto ghostStartToEndIndices(Centering const& centering,
                                               bool const includeEnd = false) const
        {
            return StartToEndIndices_(
                centering,
                [&](auto const& centering_, auto const direction) {
                    return this->ghostStartToEnd(centering_, direction);
                },
                includeEnd);
        }

        template<bool WithField, typename Centering, typename CoordsFn, typename... Indexes>
        NO_DISCARD auto inline indexesToCoords([[maybe_unused]] Centering const& centering,
                                               CoordsFn const& coordsFn,
                                               Indexes const&... indexes) const
        {
            if constexpr (WithField)
                return coordsFn(*this, centering, indexes...);
            else
                return coordsFn(*this, indexes...);
        }



        template<bool WithField = false, typename Indices, typename Centering, typename CoordsFn>
        NO_DISCARD auto indexesToCoordVectors(Indices const& indices, Centering const& centering,
                                              CoordsFn const&& coordsFn) const
        {
            auto emplace_back = [](auto& vectors, auto const&& point) {
                std::get<0>(vectors).emplace_back(point[0]);
                if constexpr (dimension > 1)
                    std::get<1>(vectors).emplace_back(point[1]);
                if constexpr (dimension > 2)
                    std::get<2>(vectors).emplace_back(point[2]);
            };

            auto xyz = tuple_fixed_type<std::vector<double>, dimension>{};

            for (auto const& indiceTuple : indices)
                std::apply(
                    [&](auto const&... args) {
                        emplace_back(
                            xyz, this->indexesToCoords<WithField>(
                                     centering, std::forward<CoordsFn const>(coordsFn), args...));
                    },
                    indiceTuple);

            return xyz;
        }


        NO_DISCARD double cellVolume() const
        {
            return std::accumulate(meshSize().begin(), meshSize().end(), 1.0,
                                   std::multiplies<double>());
        }

        /**
         * @brief physicalStartIndex returns the index of the first node of a given
         * centering and in a given direction that is in the physical domain, i.e. not a ghost node.
         */
        NO_DISCARD std::uint32_t physicalStartIndex(QtyCentering centering,
                                                    Direction direction) const
        {
            std::uint32_t icentering = static_cast<std::uint32_t>(centering);
            std::uint32_t iDir       = static_cast<std::uint32_t>(direction);
            return physicalStartIndexTable_[icentering][iDir];
        }



        NO_DISCARD std::uint32_t physicalStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                    Direction direction) const
        {
            std::uint32_t iQty                 = static_cast<std::uint32_t>(hybridQuantity);
            std::uint32_t iDir                 = static_cast<std::uint32_t>(direction);
            constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
            std::uint32_t iCentering = static_cast<std::uint32_t>(hybridQtyCentering[iQty][iDir]);

            return physicalStartIndexTable_[iCentering][iDir];
        }



        template<typename NdArrayImpl>
        NO_DISCARD std::uint32_t
        physicalStartIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                           Direction direction) const
        {
            return physicalStartIndex(field.physicalQuantity(), direction);
        }


        NO_DISCARD auto physicalStartIndex(QtyCentering centering) const
        {
            std::uint32_t icentering = static_cast<std::uint32_t>(centering);
            return physicalStartIndexTable_[icentering];
        }



        /**
         * @brief physicalEndIndex returns the index of the last node of a given
         * centering and in a given direction that is in the physical domain, i.e. not a ghost node.
         */
        NO_DISCARD std::uint32_t physicalEndIndex(QtyCentering centering, Direction direction) const
        {
            std::uint32_t icentering = static_cast<std::uint32_t>(centering);
            std::uint32_t iDir       = static_cast<std::uint32_t>(direction);

            return physicalEndIndexTable_[icentering][iDir];
        }



        NO_DISCARD std::uint32_t physicalEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                  Direction direction) const
        {
            std::uint32_t iQty                 = static_cast<std::uint32_t>(hybridQuantity);
            std::uint32_t iDir                 = static_cast<std::uint32_t>(direction);
            constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
            std::uint32_t iCentering = static_cast<std::uint32_t>(hybridQtyCentering[iQty][iDir]);

            return physicalEndIndexTable_[iCentering][iDir];
        }



        template<typename NdArrayImpl>
        NO_DISCARD std::uint32_t
        physicalEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                         Direction direction) const
        {
            return physicalEndIndex(field.physicalQuantity(), direction);
        }


        NO_DISCARD auto physicalEndIndex(QtyCentering centering) const
        {
            std::uint32_t icentering = static_cast<std::uint32_t>(centering);
            return physicalStartIndexTable_[icentering];
        }



        /**
         * @brief ghostStartIndex retuns the index of the first ghost node of a given centering
         * in a given direction. This is always zero by convention. This function exists only
         * for readibility reasons, not to have literal '0' in the code.
         */
        NO_DISCARD std::uint32_t ghostStartIndex([[maybe_unused]] QtyCentering centering,
                                                 [[maybe_unused]] Direction direction) const
        {
            // ghostStartIndex is always the first node
            return 0;
        }


        NO_DISCARD std::uint32_t
        ghostStartIndex([[maybe_unused]] HybridQuantity::Scalar const& hybridQuantity,
                        [[maybe_unused]] Direction direction) const
        {
            // ghostStartIndex is always the first node
            return 0;
        }


        template<typename NdArrayImpl>
        NO_DISCARD std::uint32_t
        ghostStartIndex([[maybe_unused]] Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                        [[maybe_unused]] Direction direction) const
        {
            // ghostStartIndex is always the first node
            return 0;
        }


        NO_DISCARD auto ghostStartIndex([[maybe_unused]] QtyCentering centering) const
        {
            return std::array<std::uint32_t, dimension>{};
        }



        /**
         * @brief ghostEndIndex returns the index of the last ghost node of a given centering
         * and in a given direction.
         */
        NO_DISCARD std::uint32_t ghostEndIndex(QtyCentering centering, Direction direction) const
        {
            std::uint32_t iCentering = static_cast<std::uint32_t>(centering);
            std::uint32_t iDir       = static_cast<std::uint32_t>(direction);

            return ghostEndIndexTable_[iCentering][iDir];
        }



        NO_DISCARD std::uint32_t ghostEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                                               Direction direction) const
        {
            std::uint32_t iQty                 = static_cast<std::uint32_t>(hybridQuantity);
            std::uint32_t iDir                 = static_cast<std::uint32_t>(direction);
            constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;
            std::uint32_t iCentering = static_cast<std::uint32_t>(hybridQtyCentering[iQty][iDir]);
            return ghostEndIndexTable_[iCentering][iDir];
        }



        template<typename NdArrayImpl>
        NO_DISCARD std::uint32_t
        ghostEndIndex(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                      Direction direction) const
        {
            return ghostEndIndex(field.physicalQuantity(), direction);
        }


        NO_DISCARD auto ghostEndIndex(QtyCentering centering) const
        {
            std::uint32_t iCentering = static_cast<std::uint32_t>(centering);
            return ghostEndIndexTable_[iCentering];
        }


        /**
         * @brief fieldNodeCoordinates returns the coordinate of a multidimensional index
         * associated with a given Field, in physical coordinates.
         */
        template<typename NdArrayImpl, typename... Indexes>
        NO_DISCARD Point<double, dimension>
        fieldNodeCoordinates(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                             const Point<double, dimension>& origin, Indexes... index) const
        {
            static_assert(sizeof...(Indexes) == dimension,
                          "Error dimension does not match number of arguments");


            std::uint32_t iQuantity       = static_cast<std::uint32_t>(field.physicalQuantity());
            constexpr std::uint32_t iDual = static_cast<std::uint32_t>(QtyCentering::dual);


            constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

            Point<std::int32_t, dimension> coord{static_cast<std::int32_t>(index)...};

            Point<double, dimension> position;

            for (std::size_t iDir = 0; iDir < dimension; ++iDir)
            {
                double halfCell = 0.0;

                auto const centering
                    = static_cast<std::uint32_t>(hybridQtyCentering[iQuantity][iDir]);
                std::int32_t const iStart = physicalStartIndexTable_[centering][iDir];

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
        NO_DISCARD Point<double, dimension> cellCenteredCoordinates(Indexes... index) const
        {
            static_assert(sizeof...(Indexes) == dimension,
                          "Error dimension does not match number of arguments");

            std::uint32_t constexpr iPrimal = static_cast<std::uint32_t>(QtyCentering::primal);

            constexpr double halfCell = 0.5;
            // A shift of +dx/2, +dy/2, +dz/2 is necessary to get the
            // cell center physical coordinates,
            // because this point is located on the dual mesh

            Point<std::uint32_t, dimension> coord(index...);

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
        NO_DISCARD auto static constexpr nbrGhosts(
            QtyCentering /*centering*/ = QtyCentering::primal)
        { // Both dual and primal ghosts are the same!
            static_assert(nbrDualGhosts_() == nbrPrimalGhosts_());

            return nbrPrimalGhosts_();
        }


        template<typename Centering, Centering centering>
        NO_DISCARD auto static constexpr nbrGhosts()
        {
            if constexpr (centering == QtyCentering::dual)
                return nbrDualGhosts_();
            else
                return nbrPrimalGhosts_();
        }



        template<typename Quantity>
        NO_DISCARD auto static constexpr nDNbrGhosts(Quantity /*centering*/ = QtyCentering::primal)
        { // Both dual and primal ghosts are the same!
            return ConstArray<std::uint32_t, dimension>(nbrGhosts());
        }


        /**
         * @brief changeCentering changes primal into dual and vice versa.
         */
        NO_DISCARD auto static changeCentering(QtyCentering centering)
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
        NO_DISCARD auto static nextIndex(QtyCentering centering, std::uint32_t indexCenter)
        {
            return indexCenter + nextIndexTable_[centering2int(centering)];
        }


        /**
         * @brief prevIndex does the same thing as nextIndex but returns the index
         * of the node of a given centering just to the left of indexCenter.
         */
        NO_DISCARD auto static prevIndex(QtyCentering centering, std::uint32_t indexCenter)
        {
            return indexCenter + prevIndexTable_[centering2int(centering)];
        }


        /** @brief returns the local 1st order derivative of the Field operand
         * at a multidimensional index and in a given direction.
         * The function can perform 1D, 2D and 3D 1st order derivatives, depending
         * on the dimensionality of the GridLayout.
         */
        template<auto direction, typename Field>
        NO_DISCARD auto deriv(Field const& operand, MeshIndex<Field::dimension> index)
        {
            auto fieldCentering = centering(operand.physicalQuantity());
            using PHARE::core::dirX;
            using PHARE::core::dirY;
            using PHARE::core::dirZ;

            if constexpr (Field::dimension == 1)
            {
                auto next = operand(nextIndex(fieldCentering[dirX], index[0]));
                auto prev = operand(prevIndex(fieldCentering[dirX], index[0]));
                return inverseMeshSize_[dirX] * (next - prev);
            }

            else if constexpr (Field::dimension == 2)
            {
                if constexpr (direction == Direction::X)
                {
                    auto next = operand(nextIndex(fieldCentering[dirX], index[0]), index[1]);
                    auto prev = operand(prevIndex(fieldCentering[dirX], index[0]), index[1]);
                    return inverseMeshSize_[dirX] * (next - prev);
                }

                if constexpr (direction == Direction::Y)
                {
                    auto next = operand(index[0], nextIndex(fieldCentering[dirY], index[1]));
                    auto prev = operand(index[0], prevIndex(fieldCentering[dirY], index[1]));
                    return inverseMeshSize_[dirY] * (next - prev);
                }
            }
            else if constexpr (Field::dimension == 3)
            {
                if constexpr (direction == Direction::X)
                {
                    auto next
                        = operand(nextIndex(fieldCentering[dirX], index[0]), index[1], index[2]);
                    auto prev
                        = operand(prevIndex(fieldCentering[dirX], index[0]), index[1], index[2]);
                    return inverseMeshSize_[dirX] * (next - prev);
                }

                if constexpr (direction == Direction::Y)
                {
                    auto next
                        = operand(index[0], nextIndex(fieldCentering[dirY], index[1]), index[2]);
                    auto prev
                        = operand(index[0], prevIndex(fieldCentering[dirY], index[1]), index[2]);
                    return inverseMeshSize_[dirY] * (next - prev);
                }
                if constexpr (direction == Direction::Z)
                {
                    auto next
                        = operand(index[0], index[1], nextIndex(fieldCentering[dirZ], index[2]));
                    auto prev
                        = operand(index[0], index[1], prevIndex(fieldCentering[dirZ], index[2]));
                    return inverseMeshSize_[dirZ] * (next - prev);
                }
            }
        }


        /** @brief returns the local laplacian of the Field operand
         * at a multidimensional index.
         * The function can perform 1D, 2D and 3D laplacian, depending
         * on the dimensionality of the GridLayout.
         */
        template<typename Field>
        NO_DISCARD auto laplacian(Field const& operand, MeshIndex<Field::dimension> index)
        {
            static_assert(Field::dimension == dimension,
                          "field dimension must be equal to gridlayout dimension");
            using PHARE::core::dirX;
            using PHARE::core::dirY;
            using PHARE::core::dirZ;

            if constexpr (Field::dimension == 1)
            {
                auto prevX = operand(index[0] - 1);
                auto hereX = operand(index[0]);
                auto nextX = operand(index[0] + 1);

                return inverseMeshSize_[dirX] * inverseMeshSize_[dirX]
                       * (nextX - 2.0 * hereX + prevX);
            }

            else if constexpr (Field::dimension == 2)
            {
                auto prevX = operand(index[0] - 1, index[1]);
                auto hereX = operand(index[0], index[1]);
                auto nextX = operand(index[0] + 1, index[1]);

                auto lapX = inverseMeshSize_[dirX] * inverseMeshSize_[dirX]
                            * (nextX - 2.0 * hereX + prevX);

                auto prevY = operand(index[0], index[1] - 1);
                auto hereY = operand(index[0], index[1]);
                auto nextY = operand(index[0], index[1] + 1);

                auto lapY = inverseMeshSize_[dirY] * inverseMeshSize_[dirY]
                            * (nextY - 2.0 * hereY + prevY);

                return lapX + lapY;
            }
            else if constexpr (Field::dimension == 3)
            {
                auto prevX = operand(index[0] - 1, index[1], index[2]);
                auto hereX = operand(index[0], index[1], index[2]);
                auto nextX = operand(index[0] + 1, index[1], index[2]);

                auto lapX = inverseMeshSize_[dirX] * inverseMeshSize_[dirX]
                            * (nextX - 2.0 * hereX + prevX);

                auto prevY = operand(index[0], index[1] - 1, index[2]);
                auto hereY = operand(index[0], index[1], index[2]);
                auto nextY = operand(index[0], index[1] + 1, index[2]);

                auto lapY = inverseMeshSize_[dirY] * inverseMeshSize_[dirY]
                            * (nextY - 2.0 * hereY + prevY);

                auto prevZ = operand(index[0], index[1], index[2] - 1);
                auto hereZ = operand(index[0], index[1], index[2]);
                auto nextZ = operand(index[0], index[1], index[2] + 1);

                auto lapZ = inverseMeshSize_[dirZ] * inverseMeshSize_[dirZ]
                            * (nextZ - 2.0 * hereZ + prevZ);

                return lapX + lapY + lapZ;
            }
        }


        /**
         * @brief localToAMR returns the AMR index associated with the given local one.
         * This method only deals with **cell** indexes.
         */
        template<typename T>
        NO_DISCARD auto localToAMR(Point<T, dimension> localPoint) const
        {
            static_assert(std::is_integral_v<T>, "Error, must be MeshIndex (integral Point)");
            Point<T, dimension> pointAMR;

            // any direction, it's the same because we want cells
            auto localStart = physicalStartIndex(QtyCentering::dual, Direction::X);

            //
            for (auto i = 0u; i < dimension; ++i)
            {
                pointAMR[i] = localPoint[i] + (AMRBox_.lower[i] - localStart);
            }
            return pointAMR;
        }


        /**
         * @brief localToAMR returns the AMR box associated with the given local one.
         * This method only deals with **cell** indexes.
         */
        template<typename T>
        NO_DISCARD auto localToAMR(Box<T, dimension> localBox) const
        {
            static_assert(std::is_integral_v<T>, "Error, must be MeshIndex (integral Point)");
            auto AMRBox = Box<T, dimension>{};

            AMRBox.lower = localToAMR(localBox.lower);
            AMRBox.upper = localToAMR(localBox.upper);

            return AMRBox;
        }


        /**
         * @brief AMRToLocal returns the local index associated with the given AMR one.
         * This method only deals with **cell** indexes.
         */
        template<typename T>
        NO_DISCARD auto AMRToLocal(Point<T, dimension> AMRPoint) const
        {
            static_assert(std::is_integral_v<T>, "Error, must be MeshIndex (integral Point)");
            Point<std::uint32_t, dimension> localPoint;

            // any direction, it's the same because we want cells
            auto localStart = physicalStartIndex(QtyCentering::dual, Direction::X);

            //
            for (auto i = 0u; i < dimension; ++i)
            {
                int local = AMRPoint[i] - (AMRBox_.lower[i] - localStart);
                assert(local >= 0);
                localPoint[i] = local;
            }
            return localPoint;
        }


        /**
         * @brief AMRToLocal returns the local Box associated with the given AMR one.
         * This method only deals with **cell** indexes.
         */
        template<typename T>
        NO_DISCARD auto AMRToLocal(Box<T, dimension> AMRBox) const
        {
            static_assert(std::is_integral_v<T>, "Error, must be MeshIndex (integral Point)");
            auto localBox = Box<std::uint32_t, dimension>{};

            localBox.lower = AMRToLocal(AMRBox.lower);
            localBox.upper = AMRToLocal(AMRBox.upper);

            return localBox;
        }



        template<typename Field, std::size_t nbr_points>
        NO_DISCARD static typename Field::type
        project(Field const& field, MeshIndex<dimension> index,
                std::array<WeightPoint<dimension>, nbr_points> wps)
        {
            typename Field::type result = 0.;
            for (auto const& wp : wps)
            {
                if constexpr (dimension == 1)
                {
                    result += wp.coef * field(index[0] + wp.indexes[0]);
                }
                if constexpr (dimension == 2)
                {
                    result += wp.coef * field(index[0] + wp.indexes[0], index[1] + wp.indexes[1]);
                }
                if constexpr (dimension == 3)
                {
                    result += wp.coef
                              * field(index[0] + wp.indexes[0], index[1] + wp.indexes[1],
                                      index[2] + wp.indexes[2]);
                }
            }
            return result;
        }



        // ----------------------------------------------------------------------
        //                      LAYOUT SPECIFIC METHODS
        //
        // The methods below return results that need information specific to the
        // layout that is used. They thus all refer to the GridLayoutImpl.
        // ----------------------------------------------------------------------


        NO_DISCARD std::string layoutName() const { return GridLayoutImpl::layoutName_; }


        /**
         * @brief returns the centering of a scalar hybrid quantity in each directions
         */
        NO_DISCARD constexpr static std::array<QtyCentering, dimension>
        centering(HybridQuantity::Scalar hybridQuantity)
        {
            return GridLayoutImpl::centering(hybridQuantity);
        }



        /**
         * @brief returns the centering of a vector hybrid quantity in each directions
         */
        NO_DISCARD constexpr static std::array<std::array<QtyCentering, dimension>, 3>
        centering(HybridQuantity::Vector hybridQuantity)
        {
            return GridLayoutImpl::centering(hybridQuantity);
        }


        /**
         * @brief GridLayout<GridLayoutImpl::dim>::allocSize
         * @return An std::array<std::uint32_t, dim> object, containing the size to which allocate
         * arrays of an HybridQuantity::Quantity 'qty' in every directions.
         */
        NO_DISCARD std::array<std::uint32_t, dimension> allocSize(HybridQuantity::Scalar qty) const
        {
            std::uint32_t iQty = static_cast<std::uint32_t>(qty);


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
        NO_DISCARD std::array<std::uint32_t, dimension> allocSizeDerived(HybridQuantity::Scalar qty,
                                                                         Direction dir) const
        {
            std::uint32_t iDerivedDir = static_cast<std::uint32_t>(dir);
            std::uint32_t iQty        = static_cast<std::uint32_t>(qty);

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
        NO_DISCARD QtyCentering
        fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction dir) const
        {
            std::uint32_t iDir = static_cast<std::uint32_t>(dir);
            std::uint32_t iQty = static_cast<std::uint32_t>(field.physicalQuantity());
            constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

            return hybridQtyCentering[iQty][iDir];
        }


        /**
         * @brief nbrPhysicalNodes returns the number of nodes in each direction, that are node
         * ghost nodes
         */
        NO_DISCARD std::array<std::uint32_t, dimension>
        nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const
        {
            std::array<QtyCentering, dimension> centerings;

            for (std::size_t iDir = 0; iDir < dimension; ++iDir)
            {
                centerings[iDir]
                    = GridLayoutImpl::hybridQtyCentering_[static_cast<std::uint32_t>(hybQty)][iDir];
            }

            return this->physicalNodeNbrFromCentering_(centerings);
        }


        /**
         * @brief derivedCentering this function returns the
         * centering (primal or dual) of a quantity after a first order derivation. dual becomes
         * primal and primal becomes dual. hybridQuantityCentering is used to know if the
         * HybridQuantity::Quantity 'qty' is primal or dual in the Direction 'dir'
         */
        NO_DISCARD QtyCentering derivedCentering(HybridQuantity::Scalar qty, Direction dir) const
        {
            std::uint32_t iField = static_cast<std::uint32_t>(qty);
            std::uint32_t idir   = static_cast<std::uint32_t>(dir);


            constexpr auto& hybridQtyCentering = GridLayoutImpl::hybridQtyCentering_;

            QtyCentering newCentering = changeCentering(hybridQtyCentering[iField][idir]);

            return newCentering;
        }


        /**
         * @brief momentsToEx return the indexes and associated coef to compute the linear
         * interpolation necessary to project moments onto Ex.
         */
        NO_DISCARD auto static constexpr momentsToEx() { return GridLayoutImpl::momentsToEx(); }


        /**
         * @brief momentsToEy return the indexes and associated coef to compute the linear
         * interpolation necessary to project moments onto Ey.
         */
        NO_DISCARD auto static constexpr momentsToEy() { return GridLayoutImpl::momentsToEy(); }


        /**
         * @brief momentsToEz return the indexes and associated coef to compute the linear
         * interpolation necessary to project moments onto Ez.
         */
        NO_DISCARD auto static constexpr momentsToEz() { return GridLayoutImpl::momentsToEz(); }



        /**
         * @brief ExToMoments return the indexes and associated coef to compute the linear
         * interpolation necessary to project Ex onto moments.
         */
        NO_DISCARD auto static constexpr ExToMoments() { return GridLayoutImpl::ExToMoments(); }



        /**
         * @brief EyToMoments return the indexes and associated coef to compute the linear
         * interpolation necessary to project Ey onto moments.
         */
        NO_DISCARD auto static constexpr EyToMoments() { return GridLayoutImpl::EyToMoments(); }



        /**
         * @brief EzToMoments return the indexes and associated coef to compute the linear
         * interpolation necessary to project Ez onto moments.
         */
        NO_DISCARD auto static constexpr EzToMoments() { return GridLayoutImpl::EzToMoments(); }


        /**
         * @brief JxToMoments return the indexes and associated coef to compute the linear
         * interpolation necessary to project Jx onto moments.
         */
        NO_DISCARD auto static constexpr JxToMoments() { return GridLayoutImpl::JxToMoments(); }


        /**
         * @brief JyToMoments return the indexes and associated coef to compute the linear
         * interpolation necessary to project Jy onto moments.
         */
        NO_DISCARD auto static constexpr JyToMoments() { return GridLayoutImpl::JyToMoments(); }


        /**
         * @brief JzToMoments return the indexes and associated coef to compute the linear
         * interpolation necessary to project Jz onto moments.
         */
        NO_DISCARD auto static constexpr JzToMoments() { return GridLayoutImpl::JzToMoments(); }


        /**
         * @brief ByToEx return the indexes and associated coef to compute the linear
         * interpolation necessary to project By onto Ex.
         */
        NO_DISCARD auto static constexpr ByToEx() { return GridLayoutImpl::ByToEx(); }


        /**
         * @brief BzToEx return the indexes and associated coef to compute the linear
         * interpolation necessary to project Bz onto Ex.
         */
        NO_DISCARD auto static constexpr BzToEx() { return GridLayoutImpl::BzToEx(); }


        /**
         * @brief BxToEy return the indexes and associated coef to compute the linear
         * interpolation necessary to project Bx onto Ey.
         */
        NO_DISCARD auto static constexpr BxToEy() { return GridLayoutImpl::BxToEy(); }



        /**
         * @brief BzToEy return the indexes and associated coef to compute the linear
         * interpolation necessary to project Bz onto Ey.
         */
        NO_DISCARD auto static constexpr BzToEy() { return GridLayoutImpl::BzToEy(); }



        /**
         * @brief BxToEz return the indexes and associated coef to compute the linear
         * interpolation necessary to project Bx onto Ez.
         */
        NO_DISCARD auto static constexpr BxToEz() { return GridLayoutImpl::BxToEz(); }



        /**
         * @brief ByToEz return the indexes and associated coef to compute the linear
         * interpolation necessary to project By onto Ez.
         */
        NO_DISCARD auto static constexpr ByToEz() { return GridLayoutImpl::ByToEz(); }



        /**
         * @brief JxToEx return the indexes and associated coef to compute the linear
         * interpolation necessary to project Jx onto Ex.
         */
        NO_DISCARD auto static constexpr JxToEx() { return GridLayoutImpl::JxToEx(); }



        /**
         * @brief JyToEy return the indexes and associated coef to compute the linear
         * interpolation necessary to project Jy onto Ey.
         */
        NO_DISCARD auto static constexpr JyToEy() { return GridLayoutImpl::JyToEy(); }



        /**
         * @brief JzToEz return the indexes and associated coef to compute the linear
         * interpolation necessary to project Jz onto Ez.
         */
        NO_DISCARD auto static constexpr JzToEz() { return GridLayoutImpl::JzToEz(); }




        template<typename Field, typename Fn>
        void evalOnBox(Field& field, Fn&& fn) const
        {
            auto indices = [&](auto const& centering, auto const direction) {
                return this->physicalStartToEnd(centering, direction);
            };

            evalOnBox_(field, fn, indices);
        }

        template<typename Field, typename Fn>
        void evalOnGhostBox(Field& field, Fn&& fn) const
        {
            auto indices = [&](auto const& centering, auto const direction) {
                return this->ghostStartToEnd(centering, direction);
            };

            evalOnBox_(field, fn, indices);
        }


    private:
        template<typename Field, typename IndicesFn, typename Fn>
        static void evalOnBox_(Field& field, Fn& fn, IndicesFn& startToEnd)
        {
            auto const [ix0, ix1] = startToEnd(field, Direction::X);
            for (auto ix = ix0; ix <= ix1; ++ix)
            {
                if constexpr (dimension == 1)
                {
                    fn(ix);
                }
                else
                {
                    auto const [iy0, iy1] = startToEnd(field, Direction::Y);

                    for (auto iy = iy0; iy <= iy1; ++iy)
                    {
                        if constexpr (dimension == 2)
                        {
                            fn(ix, iy);
                        }
                        else
                        {
                            auto const [iz0, iz1] = startToEnd(field, Direction::Z);

                            for (auto iz = iz0; iz <= iz1; ++iz)
                                fn(ix, iy, iz);
                        }
                    }
                }
            }
        }


        template<typename Centering, typename StartToEnd>
        auto StartToEndIndices_(Centering const& centering, StartToEnd const&& startToEnd,
                                bool const includeEnd = false) const
        {
            std::vector<tuple_fixed_type<std::uint32_t, dimension>> indices;

            std::size_t plus = (includeEnd) ? 1 : 0;

            auto const [ix0, ix1] = startToEnd(centering, Direction::X);

            for (auto ix = ix0; ix < ix1 + plus; ++ix)
            {
                if constexpr (dimension > 1)
                {
                    auto const [iy0, iy1] = startToEnd(centering, Direction::Y);

                    for (auto iy = iy0; iy < iy1 + plus; ++iy)
                    {
                        if constexpr (dimension > 2)
                        {
                            auto const [iz0, iz1] = startToEnd(centering, Direction::Z);

                            for (auto iz = iz0; iz < iz1 + plus; ++iz)
                                indices.emplace_back(std::make_tuple(ix, iy, iz));
                        }
                        else // 2D
                        {
                            indices.emplace_back(std::make_tuple(ix, iy));
                        }
                    }
                }
                else // 1D
                {
                    indices.emplace_back(std::make_tuple(ix));
                }
            }
            return indices;
        }




        /**
         * @brief nextPrimal_ returns the index shift needed to go to the next primal
         * node from a dual node. This depends on whether the dual have more ghost nodes
         * than the primal.
         * returning 0 means that the next primal has the same index as the current dual.
         * returning 1 means that the next primal index is the current dual + 1
         */
        NO_DISCARD constexpr static auto nextPrimal_()
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
        NO_DISCARD constexpr static auto prevPrimal_()
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
        NO_DISCARD constexpr static auto nextDual_()
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
        NO_DISCARD constexpr static auto prevDual_()
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
        NO_DISCARD std::uint32_t constexpr static nbrDualGhosts_()
        {
            return (interp_order + 1) / 2 + nbrParticleGhosts();
        }


        /**
         * @brief nbrPrimalGhosts_ returns the number of primal ghost nodes.
         * Contrary to dual ghost nodes, the formula to get the number of primal ghost nodes depend
         * on the interpolation order. If based only on the particle-mesh interaction, order1 would
         * not need primal ghost nodes. But there is a minimum of 1 ghost node for primal that is
         * linked to the possibility of calculating second order derivatives of primal quantities
         * (e.g. laplacian of J for a yee lattice). Dual ghosts don't have this issue since they
         * always have at least 1 ghost.
         */
        NO_DISCARD std::uint32_t constexpr static nbrPrimalGhosts_() { return nbrDualGhosts_(); }



        NO_DISCARD std::uint32_t static constexpr dualOffset_() noexcept { return 1; }


        /**
         * @brief physicalNodeNbrFromCentering_ returns the number of physical nodes for all
         * directions depending on the multi-dimensional centering.
         */
        NO_DISCARD std::array<std::uint32_t, dimension> physicalNodeNbrFromCentering_(
            std::array<QtyCentering, dimension> const& qtyCenterings) const
        {
            std::array<std::uint32_t, dimension> nodeNbr;

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
        NO_DISCARD std::array<std::uint32_t, dimension>
        nodeNbrFromCentering_(std::array<QtyCentering, dimension> const& qtyCenterings) const
        {
            std::array<std::uint32_t, dimension> nbrNodes
                = physicalNodeNbrFromCentering_(qtyCenterings);

            for (std::size_t iDir = 0; iDir < dimension; ++iDir)
            {
                nbrNodes[iDir] += 2 * nbrGhosts(qtyCenterings[iDir]);
            }

            return nbrNodes;
        }


        NO_DISCARD auto initPhysicalStart_()
        {
            std::array<std::array<std::uint32_t, dimension>, 2> physicalStartIndexTable;

            std::uint32_t iprimal = static_cast<std::uint32_t>(data.primal);
            std::uint32_t idual   = static_cast<std::uint32_t>(data.dual);

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
        NO_DISCARD auto initPhysicalEnd_()
        {
            std::array<std::array<std::uint32_t, dimension>, 2> physicalEndIndexTable;

            std::uint32_t iprimal = static_cast<std::uint32_t>(data.primal);
            std::uint32_t idual   = static_cast<std::uint32_t>(data.dual);

            physicalEndIndexTable[iprimal][data.idirX]
                = physicalStartIndexTable_[iprimal][data.idirX] + nbrPhysicalCells_[data.idirX];


            physicalEndIndexTable[idual][data.idirX] = physicalStartIndexTable_[idual][data.idirX]
                                                       + nbrPhysicalCells_[data.idirX]
                                                       - dualOffset_();

            if constexpr (dimension > 1)
            {
                physicalEndIndexTable[iprimal][data.idirY]
                    = physicalStartIndexTable_[iprimal][data.idirY] + nbrPhysicalCells_[data.idirY];

                physicalEndIndexTable[idual][data.idirY]
                    = physicalStartIndexTable_[idual][data.idirY] + nbrPhysicalCells_[data.idirY]
                      - dualOffset_();

                if constexpr (dimension > 2)
                {
                    physicalEndIndexTable[iprimal][data.idirZ]
                        = physicalStartIndexTable_[iprimal][data.idirZ]
                          + nbrPhysicalCells_[data.idirZ];

                    physicalEndIndexTable[idual][data.idirZ]
                        = physicalStartIndexTable_[idual][data.idirZ]
                          + nbrPhysicalCells_[data.idirZ] - dualOffset_();
                }
            }
            return physicalEndIndexTable;
        }



        /**
         * @brief GridLayout<GridLayoutImpl::dim>::initGhostEnd calculate and stores the index
         * of the last primal and dual nodes in each direction. The formula simply
         * consists in starting at physicalEndIndex() and to add the number of ghost nodes.
         */
        NO_DISCARD auto initGhostEnd_()
        {
            std::array<std::array<std::uint32_t, dimension>, 2> ghostEndIndexTable;

            std::uint32_t iprimal = static_cast<std::uint32_t>(data.primal);
            std::uint32_t idual   = static_cast<std::uint32_t>(data.dual);

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
        std::array<std::uint32_t, dimension> nbrPhysicalCells_;
        std::array<double, dimension> inverseMeshSize_;
        static constexpr gridDataT data{};

        // stores key indices in each direction (3) for primal and dual nodes (2)
        std::array<std::array<std::uint32_t, dimension>, 2> physicalStartIndexTable_;
        std::array<std::array<std::uint32_t, dimension>, 2> physicalEndIndexTable_;
        std::array<std::array<std::uint32_t, dimension>, 2> ghostEndIndexTable_;
        Box<int, dimension> AMRBox_;


        // this constexpr initialization only works if primal==0 and dual==1
        // this is defined in gridlayoutdefs.hpp don't change it because these
        // arrays will be accessed with [primal] and [dual] indexes.
        constexpr static std::array<int, 2> nextIndexTable_{{nextPrimal_(), nextDual_()}};
        constexpr static std::array<int, 2> prevIndexTable_{{prevPrimal_(), prevDual_()}};
    };


} // namespace core
} // namespace PHARE

#endif // GridLayout_HPP
