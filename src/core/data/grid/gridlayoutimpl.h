#ifndef PHARE_CORE_GRID_GRIDLAYOUTIMPL_H
#define PHARE_CORE_GRID_GRIDLAYOUTIMPL_H



#include "data/field/field.h"
#include "hybrid/hybrid_quantities.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"
#include "utilities/types.h"

#include "gridlayoutdefs.h"

#include <array>
#include <vector>

namespace PHARE
{
struct WeightPoint
{
    int32 ix, iy, iz;
    double coef;
};


using LinearCombination = std::vector<WeightPoint>;

enum class Layout { Yee };


/**
 * @brief GridLayoutImpl is a virtual pure class, it declares interface
 * prototypes publically inherited by GridLayoutImplYee.
 *
 * The virtual methods are related to grid layout operations:
 * - physical domain start/end indexes
 * - indexes of the first and last ghost nodes
 * - allocation sizes for Field attributes of other classes
 * - partial derivative operator (Faraday)
 * - physical coordinate given a field and a primal point (ix, iy, iz)
 * - cell centered coordinate given a primal point (ix, iy, iz)
 *
 */
template<Layout, std::size_t dim>
class GridLayoutImpl
{
public:
    constexpr static std::array<QtyCentering, dim>
    centering(HybridQuantity::Scalar const& hybridQuantity);

    // start and end index used in computing loops
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
    void deriv1D(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
                 Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const;
    /* template<typename NdArrayImpl> */
    /*   void deriv2D(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand, */
    /*                Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const; */
    /* template<typename NdArrayImpl> */
    /*   void deriv3D(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand, */
    /*                Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const; */

    template<typename NdArrayImpl, typename... Indexes>
    Point<double, dim> fieldNodeCoordinates(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                                            const Point<double, dim>& origin,
                                            Indexes... index) const;

    template<typename... Indexes>
    Point<double, dim> cellCenteredCoordinates(Indexes... index) const;

    template<typename NdArrayImpl>
    QtyCentering fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                Direction dir);

    uint32 nbrGhostNodes(QtyCentering centering) const;

    template<typename NdArrayImpl>
    std::array<uint32, dim>
    nbrPhysicalNodes(Field<NdArrayImpl, HybridQuantity::Scalar> const& field);

    std::array<uint32, dim> nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const;
    uint32 nbDimensions() const;



    LinearCombination const& momentsToEx() const;
    LinearCombination const& momentsToEy() const;
    LinearCombination const& momentsToEz() const;

    LinearCombination const& ByToEx() const;
    LinearCombination const& ByToEz() const;

    LinearCombination const& BxToEy() const;
    LinearCombination const& BxToEz() const;

    LinearCombination const& BzToEx() const;
    LinearCombination const& BzToEy() const;

    LinearCombination const& ExToMoment() const;
    LinearCombination const& EyToMoment() const;
    LinearCombination const& EzToMoment() const;
};
} // namespace PHARE

#endif // GRIDLAYOUTIMPL_H
