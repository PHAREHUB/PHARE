#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_H
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_H



#include "gridlayoutdefs.h"
#include "gridlayoutimpl.h"
#include "gridlayoutimplinternals.h"

#include "hybrid/hybrid_quantities.h"
#include "utilities/constants.h"
#include "utilities/types.h"

#include <array>
#include <vector>

namespace PHARE
{
/**
 * @brief GridLayoutNdArrayImplYee class is a concrete GridLayoutNdArrayImpl used a Yee
 * type grid layout is needed.
 *
 * It provides methods related to grid layout operations:
 * - physical domain start/end indexes
 * - indexes of the first and last ghost nodes
 * - allocation sizes for Field attributes of other classes
 * - partial derivative operator (Faraday)
 * - physical coordinate given a field and a primal point (ix, iy, iz)
 * - cell centered coordinate given a primal point (ix, iy, iz)
 */
template<std::size_t dim>
class GridLayoutImpl<Layout::Yee, dim>
    : private GridLayoutImplInternals<GridLayoutImpl<Layout::Yee, dim>, dim>
{
    friend class GridLayoutImplInternals<GridLayoutImpl<Layout::Yee, dim>, dim>;

    // ------------------------------------------------------------------------
    //                              PRIVATE
    // ------------------------------------------------------------------------
private:
    void initLinearCombinations_();

    LinearCombination momentsToEx_;
    LinearCombination momentsToEy_;
    LinearCombination momentsToEz_;
    LinearCombination BxToEy_;
    LinearCombination BxToEz_;
    LinearCombination ByToEx_;
    LinearCombination ByToEz_;
    LinearCombination BzToEx_;
    LinearCombination BzToEy_;
    LinearCombination ExToMoment_;
    LinearCombination EyToMoment_;
    LinearCombination EzToMoment_;

    /**
     * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::initLayoutCentering_ initialize the
     * table hybridQuantityCentering_. This is THE important array in the GridLayout module. This
     * table knows which quantity is primal/dual along each direction. It is **this** array that
     * **defines** what a Yee Layout is. Once this array is defined, the rest of the GridLayout
     * needs this array OK and can go on from here... hence all other functions in the Yee interface
     * are just calling private implementation common to all layouts
     */
    constexpr auto static initLayoutCentering_()
    {
        const gridDataT data{};
        const std::array<QtyCentering, NBR_COMPO> Bx = {{data.primal, data.dual, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> By = {{data.dual, data.primal, data.dual}};
        const std::array<QtyCentering, NBR_COMPO> Bz = {{data.dual, data.dual, data.primal}};

        const std::array<QtyCentering, NBR_COMPO> Ex = {{data.dual, data.primal, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Ey = {{data.primal, data.dual, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Ez = {{data.primal, data.primal, data.dual}};


        const std::array<QtyCentering, NBR_COMPO> Jx = {{data.dual, data.primal, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Jy = {{data.primal, data.dual, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Jz = {{data.primal, data.primal, data.dual}};

        const std::array<QtyCentering, NBR_COMPO> Rho = {{data.primal, data.primal, data.primal}};

        const std::array<QtyCentering, NBR_COMPO> Vx = {{data.primal, data.primal, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Vy = {{data.primal, data.primal, data.primal}};
        const std::array<QtyCentering, NBR_COMPO> Vz = {{data.primal, data.primal, data.primal}};

        const std::array<QtyCentering, NBR_COMPO> P = {{data.primal, data.primal, data.primal}};

        const std::array<std::array<QtyCentering, NBR_COMPO>,
                         static_cast<std::size_t>(HybridQuantity::Scalar::count)>
            hybridQtyCentering{Bx, By, Bz, Ex, Ey, Ez, Jx, Jy, Jz, Rho, Vx, Vy, Vz, P};


        return hybridQtyCentering;
    }
    //! says for each HybridQuantity::Quantity whether it is primal or dual, in each direction
    constexpr const static std::array<std::array<QtyCentering, NBR_COMPO>,
                                      static_cast<std::size_t>(HybridQuantity::Scalar::count)>
        hybridQtyCentering_{initLayoutCentering_()};

    static const std::size_t dim_{dim};

    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------
public:
    GridLayoutImpl(std::array<double, dim> const& dx, std::array<uint32, dim> const& nbrCells,
                   Point<double, dim> const& origin, uint32 interpOrder);


    ~GridLayoutImpl() = default;

    constexpr static uint32 nbDimensions() { return dim; }


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



    template<typename NdArrayImpl>
    void deriv1D(Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
                 Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const;

    std::array<uint32, dim> allocSize(HybridQuantity::Scalar qtyType) const;

    std::array<uint32, dim> allocSizeDerived(HybridQuantity::Scalar qty, Direction dir) const;

    template<typename NdArrayImpl, typename... Indexes>
    Point<double, dim> fieldNodeCoordinates(const Field<NdArrayImpl, HybridQuantity::Scalar>& field,
                                            const Point<double, dim>& origin,
                                            Indexes... index) const;

    template<typename... Indexes>
    Point<double, dim> cellCenteredCoordinates(Indexes... index) const;

    template<typename NdArrayImpl>
    QtyCentering fieldCentering(Field<NdArrayImpl, HybridQuantity::Scalar> const& field,
                                Direction dir) const;

    uint32 nbrGhostNodes(QtyCentering centering) const;

    template<typename NdArrayImpl>
    std::array<uint32, dim>
    nbrPhysicalNodes(Field<NdArrayImpl, HybridQuantity::Scalar> const& field);

    std::array<uint32, dim> nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const;

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




template<std::size_t dim>
GridLayoutImpl<Layout::Yee, dim>::GridLayoutImpl(std::array<double, dim> const& dx,
                                                 std::array<uint32, dim> const& nbrCells,
                                                 Point<double, dim> const& origin,
                                                 uint32 interpOrder)
    : GridLayoutImplInternals<GridLayoutImpl<Layout::Yee, dim>, dim>(dx, nbrCells, interpOrder,
                                                                     origin)
{
    constexpr gridDataT gridData{};
    this->initPhysicalStart(gridData);

    this->initPhysicalEnd(gridData);
    this->initGhostEnd(gridData);
    this->initLinearCombinations_();
}



template<>
constexpr std::array<QtyCentering, 1>
GridLayoutImpl<Layout::Yee, 1>::centering(HybridQuantity::Scalar const& hybridQuantity)
{
    constexpr gridDataT gridData_{};
    switch (hybridQuantity)
    {
        case HybridQuantity::Scalar::Bx:
            return {{hybridQtyCentering_[gridData_.iBx][gridData_.idirX]}};
        case HybridQuantity::Scalar::By:
            return {{hybridQtyCentering_[gridData_.iBy][gridData_.idirX]}};
        case HybridQuantity::Scalar::Bz:
            return {{hybridQtyCentering_[gridData_.iBz][gridData_.idirX]}};
        case HybridQuantity::Scalar::Ex:
            return {{hybridQtyCentering_[gridData_.iEx][gridData_.idirX]}};
        case HybridQuantity::Scalar::Ey:
            return {{hybridQtyCentering_[gridData_.iEy][gridData_.idirX]}};
        case HybridQuantity::Scalar::Ez:
            return {{hybridQtyCentering_[gridData_.iEz][gridData_.idirX]}};
        case HybridQuantity::Scalar::Jx:
            return {{hybridQtyCentering_[gridData_.iJx][gridData_.idirX]}};
        case HybridQuantity::Scalar::Jy:
            return {{hybridQtyCentering_[gridData_.iJy][gridData_.idirX]}};
        case HybridQuantity::Scalar::Jz:
            return {{hybridQtyCentering_[gridData_.iJz][gridData_.idirX]}};
        case HybridQuantity::Scalar::rho:
            return {{hybridQtyCentering_[gridData_.irho][gridData_.idirX]}};
        case HybridQuantity::Scalar::Vx:
            return {{hybridQtyCentering_[gridData_.iVx][gridData_.idirX]}};
        case HybridQuantity::Scalar::Vy:
            return {{hybridQtyCentering_[gridData_.iVy][gridData_.idirX]}};
        case HybridQuantity::Scalar::Vz:
            return {{hybridQtyCentering_[gridData_.iVz][gridData_.idirX]}};
        case HybridQuantity::Scalar::P:
            return {{hybridQtyCentering_[gridData_.iP][gridData_.idirX]}};
        default: throw std::runtime_error("Wrong hybridQuantity");
    }
}
template<>
constexpr std::array<QtyCentering, 2>
GridLayoutImpl<Layout::Yee, 2>::centering(HybridQuantity::Scalar const& hybridQuantity)
{
    constexpr gridDataT gridData_{};
    switch (hybridQuantity)
    {
        case HybridQuantity::Scalar::Bx:
            return {{hybridQtyCentering_[gridData_.iBx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iBx][gridData_.idirY]}};
        case HybridQuantity::Scalar::By:
            return {{hybridQtyCentering_[gridData_.iBy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iBy][gridData_.idirY]}};
        case HybridQuantity::Scalar::Bz:
            return {{hybridQtyCentering_[gridData_.iBz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iBz][gridData_.idirY]}};
        case HybridQuantity::Scalar::Ex:
            return {{hybridQtyCentering_[gridData_.iEx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iEx][gridData_.idirY]}};
        case HybridQuantity::Scalar::Ey:
            return {{hybridQtyCentering_[gridData_.iEy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iEy][gridData_.idirY]}};
        case HybridQuantity::Scalar::Ez:
            return {{hybridQtyCentering_[gridData_.iEz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iEz][gridData_.idirY]}};
        case HybridQuantity::Scalar::Jx:
            return {{hybridQtyCentering_[gridData_.iJx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iJx][gridData_.idirY]}};
        case HybridQuantity::Scalar::Jy:
            return {{hybridQtyCentering_[gridData_.iJy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iJy][gridData_.idirY]}};
        case HybridQuantity::Scalar::Jz:
            return {{hybridQtyCentering_[gridData_.iJz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iJz][gridData_.idirY]}};
        case HybridQuantity::Scalar::rho:
            return {{hybridQtyCentering_[gridData_.irho][gridData_.idirX],
                     hybridQtyCentering_[gridData_.irho][gridData_.idirY]}};
        case HybridQuantity::Scalar::Vx:
            return {{hybridQtyCentering_[gridData_.iVx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iVx][gridData_.idirY]}};
        case HybridQuantity::Scalar::Vy:
            return {{hybridQtyCentering_[gridData_.iVy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iVy][gridData_.idirY]}};
        case HybridQuantity::Scalar::Vz:
            return {{hybridQtyCentering_[gridData_.iVz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iVz][gridData_.idirY]}};
        case HybridQuantity::Scalar::P:
            return {{hybridQtyCentering_[gridData_.iP][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iP][gridData_.idirY]}};
        default: throw std::runtime_error("Wrong hybridQuantity");
    }
}




template<>
constexpr std::array<QtyCentering, 3>
GridLayoutImpl<Layout::Yee, 3>::centering(HybridQuantity::Scalar const& hybridQuantity)
{
    constexpr gridDataT gridData_{};
    switch (hybridQuantity)
    {
        case HybridQuantity::Scalar::Bx:
            return {{hybridQtyCentering_[gridData_.iBx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iBx][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iBx][gridData_.idirZ]}};
        case HybridQuantity::Scalar::By:
            return {{hybridQtyCentering_[gridData_.iBy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iBy][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iBy][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Bz:
            return {{hybridQtyCentering_[gridData_.iBz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iBz][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iBz][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Ex:
            return {{hybridQtyCentering_[gridData_.iEx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iEx][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iEx][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Ey:
            return {{hybridQtyCentering_[gridData_.iEy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iEy][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iEy][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Ez:
            return {{hybridQtyCentering_[gridData_.iEz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iEz][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iEz][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Jx:
            return {{hybridQtyCentering_[gridData_.iJx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iJx][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iJx][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Jy:
            return {{hybridQtyCentering_[gridData_.iJy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iJy][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iJy][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Jz:
            return {{hybridQtyCentering_[gridData_.iJz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iJz][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iJz][gridData_.idirZ]}};
        case HybridQuantity::Scalar::rho:
            return {{hybridQtyCentering_[gridData_.irho][gridData_.idirX],
                     hybridQtyCentering_[gridData_.irho][gridData_.idirY],
                     hybridQtyCentering_[gridData_.irho][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Vx:
            return {{hybridQtyCentering_[gridData_.iVx][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iVx][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iVx][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Vy:
            return {{hybridQtyCentering_[gridData_.iVy][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iVy][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iVy][gridData_.idirZ]}};
        case HybridQuantity::Scalar::Vz:
            return {{hybridQtyCentering_[gridData_.iVz][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iVz][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iVz][gridData_.idirZ]}};
        case HybridQuantity::Scalar::P:
            return {{hybridQtyCentering_[gridData_.iP][gridData_.idirX],
                     hybridQtyCentering_[gridData_.iP][gridData_.idirY],
                     hybridQtyCentering_[gridData_.iP][gridData_.idirZ]}};
        default: throw std::runtime_error("Wrong hybridQuantity");
    }
}




template<std::size_t dim>
void GridLayoutImpl<Layout::Yee, dim>::initLinearCombinations_()
{
    // cf https://hephaistos.lpp.polytechnique.fr/redmine/projects/hyb-par/wiki/Ohm
    // for how to calculate coefficients and shift indexes.

    int dualToPrimal;
    int primalTodual;

    if (this->interpOrder_ == 1 || this->interpOrder_ == 2 || this->interpOrder_ == 4)
    {
        dualToPrimal = -1;
        primalTodual = 1;
    }
    else if (this->interpOrder_ == 3)
    {
        dualToPrimal = 1;
        primalTodual = -1;
    }
    else
    {
        throw std::runtime_error("GridLayout Yee cannot be initialized: wrong interpolation order");
    }


    WeightPoint P1;
    WeightPoint P2;

    // moment to Ex is Ppp to Dpp
    // shift only in X
    // the average is done for all simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = 0.5;
    P2.ix   = primalTodual;
    P2.iy   = 0;
    P2.iz   = 0;
    P2.coef = 0.5;
    momentsToEx_.push_back(P1);
    momentsToEx_.push_back(P2);


    // moment to Ey is pPp to pDp
    // shift only in Y
    // the average is done only for 2D and 3D simulation
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim >= 2) ? 0.5 : 1.;
    momentsToEy_.push_back(P1);

    // in 2 and 3D, add another point and average
    if (dim >= 2)
    {
        P2.ix   = 0;
        P2.iy   = primalTodual;
        P2.iz   = 0;
        P2.coef = 0.5;
        momentsToEy_.push_back(P2);
    }


    // moment to Ez is ppP to ppD
    // shift only in Z
    // the average is done only for 3D simulation
    // hence for 1D and 2D runs coef==1
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim == 3) ? 0.5 : 1;
    momentsToEz_.push_back(P1);

    if (dim == 3)
    {
        P2.ix   = 0;
        P2.iy   = 0;
        P2.iz   = primalTodual;
        P2.coef = 0.5;
        momentsToEz_.push_back(P2);
    }




    // Bx to Ey is pdD to pdP
    // shift only in Z
    // the average is done only for 3D simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim == 3) ? 0.5 : 1;
    BxToEy_.push_back(P1);

    if (dim == 3)
    {
        P2.ix   = 0;
        P2.iy   = 0;
        P2.iz   = dualToPrimal;
        P2.coef = 0.5;
        BxToEy_.push_back(P2);
    }



    // Bx to Ez is pDd to pPd
    // shift in the Y direction only
    // the average is done for 2D and 3D simulations
    // hence for 1D simulations coef is 1
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim >= 2) ? 0.5 : 1;
    BxToEz_.push_back(P1);

    if (dim >= 2)
    {
        P2.ix   = 0;
        P2.iy   = dualToPrimal;
        P2.iz   = 0;
        P2.coef = 0.5;
        BxToEz_.push_back(P2);
    }


    // By to Ex is dpD to dpP
    // shift only in the Z direction
    // averaging is done only for 3D simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim == 3) ? 0.5 : 1;
    ByToEx_.push_back(P1);

    if (dim == 3)
    {
        P2.ix   = 0;
        P2.iy   = 0;
        P2.iz   = dualToPrimal;
        P2.coef = 0.5;
        ByToEx_.push_back(P2);
    }

    // By to Ez is Dpd to Ppd
    // shift only in the X direction
    // the averaging is done in all simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = 0.5;
    P2.ix   = dualToPrimal;
    P2.iy   = 0;
    P2.iz   = 0;
    P2.coef = 0.5;
    ByToEz_.push_back(P1);
    ByToEz_.push_back(P2);


    // Bz to Ex is dDp to dPp
    // shift only in the Y direction
    // the averaging is done for 2D and 3D simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim >= 2) ? 0.5 : 1;
    BzToEx_.push_back(P1);

    if (dim >= 2)
    {
        P2.ix   = 0;
        P2.iy   = dualToPrimal;
        P2.iz   = 0;
        P2.coef = 0.5;
        BzToEx_.push_back(P2);
    }


    // Bz to Ey is Ddp to Pdp
    // shift only in the X direction
    // the averaging is done for all simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = 0.5;
    BzToEy_.push_back(P1);
    P2.ix   = dualToPrimal;
    P2.iy   = 0;
    P2.iz   = 0;
    P2.coef = 0.5;
    BzToEy_.push_back(P2);


    // Ex to Moment is Dpp to Ppp
    // shift only in the X direction
    // the averaging is done for all simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = 0.5;
    ExToMoment_.push_back(P1);
    P2.ix   = dualToPrimal;
    P2.iy   = 0;
    P2.iz   = 0;
    P2.coef = 0.5;
    ExToMoment_.push_back(P2);


    // Ey to Moment is pDp to PPP
    // shift tis only in the Y direction
    // the averaging is done for 2D and 3D simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim >= 2) ? 0.5 : 1;
    EyToMoment_.push_back(P1);

    if (dim >= 2)
    {
        P2.ix   = 0;
        P2.iy   = dualToPrimal;
        P2.iz   = 0;
        P2.coef = 0.5;
        EyToMoment_.push_back(P2);
    }


    // Ez to Moment is ppD on ppP
    // shift only in the Z direction
    // the averaging is only for 3D simulations
    P1.ix   = 0;
    P1.iy   = 0;
    P1.iz   = 0;
    P1.coef = (dim == 3) ? 0.5 : 1;
    EzToMoment_.push_back(P1);

    if (dim == 3)
    {
        P2.ix   = 0;
        P2.iy   = 0;
        P2.iz   = dualToPrimal;
        P2.coef = 0.5;
        EzToMoment_.push_back(P2);
    }
}




template<std::size_t dim>
std::array<uint32, dim>
GridLayoutImpl<Layout::Yee, dim>::allocSize(HybridQuantity::Scalar qty) const
{
    return this->allocSize_(qty);
}




// TODO : WARNING 1st order only
// Can it be moved to ImplInternals?
template<std::size_t dim>
std::array<uint32, dim>
GridLayoutImpl<Layout::Yee, dim>::allocSizeDerived(HybridQuantity::Scalar qty, Direction dir) const
{
    return this->allocSizeDerived_(qty, dir);
}




// start and end index used in computing loops
template<std::size_t dim>
uint32 GridLayoutImpl<Layout::Yee, dim>::physicalStartIndex(QtyCentering centering,
                                                            Direction direction) const
{
    return this->physicalStartIndex_(centering, direction);
}


template<std::size_t dim>
uint32
GridLayoutImpl<Layout::Yee, dim>::physicalStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                     Direction direction) const
{
    return this->physicalStartIndex_(hybridQuantity, direction);
}

template<std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImpl<Layout::Yee, dim>::physicalStartIndex(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    return this->physicalStartIndex_(field, direction);
}


template<std::size_t dim>
uint32 GridLayoutImpl<Layout::Yee, dim>::physicalEndIndex(QtyCentering centering,
                                                          Direction direction) const
{
    return this->physicalEndIndex_(centering, direction);
}


template<std::size_t dim>
uint32
GridLayoutImpl<Layout::Yee, dim>::physicalEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                   Direction direction) const
{
    return this->physicalEndIndex_(hybridQuantity, direction);
}


template<std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImpl<Layout::Yee, dim>::physicalEndIndex(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    return this->physicalEndIndex_(field, direction);
}


template<std::size_t dim>
uint32 GridLayoutImpl<Layout::Yee, dim>::ghostStartIndex(QtyCentering centering,
                                                         Direction direction) const
{
    // should we directly return 0 and remove ghostStartIndex_ ?
    return this->ghostStartIndex_(centering, direction);
}

template<std::size_t dim>
uint32
GridLayoutImpl<Layout::Yee, dim>::ghostStartIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                  Direction direction) const
{
    // should we directly return 0 and remove ghostStartIndex_ ?
    return this->ghostStartIndex_(hybridQuantity, direction);
}

template<std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImpl<Layout::Yee, dim>::ghostStartIndex(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    // should we directly return 0 and remove ghostStartIndex_ ?
    return this->ghostStartIndex_(field, direction);
}



template<std::size_t dim>
uint32 GridLayoutImpl<Layout::Yee, dim>::ghostEndIndex(QtyCentering centering,
                                                       Direction direction) const
{
    return this->ghostEndIndex_(centering, direction);
}

template<std::size_t dim>
uint32 GridLayoutImpl<Layout::Yee, dim>::ghostEndIndex(HybridQuantity::Scalar const& hybridQuantity,
                                                       Direction direction) const
{
    return this->ghostEndIndex_(hybridQuantity, direction);
}



template<std::size_t dim>
template<typename NdArrayImpl>
uint32 GridLayoutImpl<Layout::Yee, dim>::ghostEndIndex(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction direction) const
{
    return this->ghostEndIndex_(field, direction);
}


template<std::size_t dim>
template<typename NdArrayImpl, typename... Indexes>
Point<double, dim> GridLayoutImpl<Layout::Yee, dim>::fieldNodeCoordinates(
    const Field<NdArrayImpl, HybridQuantity::Scalar>& field, const Point<double, dim>& origin,
    Indexes... index) const
{
    return this->fieldNodeCoordinates_(field, origin, index...);
}




template<std::size_t dim>
template<typename... Indexes>
Point<double, dim> GridLayoutImpl<Layout::Yee, dim>::cellCenteredCoordinates(Indexes... index) const
{
    return this->cellCenteredCoordinates_(index...);
}



template<std::size_t dim>
template<typename NdArrayImpl>
QtyCentering GridLayoutImpl<Layout::Yee, dim>::fieldCentering(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& field, Direction dir) const
{
    return this->fieldCentering_(field, dir);
}


template<std::size_t dim>
uint32 GridLayoutImpl<Layout::Yee, dim>::nbrGhostNodes(QtyCentering centering) const
{
    return this->nbrGhosts(centering);
}




template<std::size_t dim>
std::array<uint32, dim>
GridLayoutImpl<Layout::Yee, dim>::nbrPhysicalNodes(HybridQuantity::Scalar hybQty) const
{
    std::array<QtyCentering, dim> centerings;

    for (std::size_t iDir = 0; iDir < dim; ++iDir)
    {
        centerings[iDir] = this->hybridQtyCentering_[static_cast<uint32>(hybQty)][iDir];
    }

    return this->physicalNodeNbrFromCentering_(centerings);
}




/**
 * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::deriv1D It was decided to compute the
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
template<std::size_t dim>
template<typename NdArrayImpl>
void GridLayoutImpl<Layout::Yee, dim>::deriv1D(
    Field<NdArrayImpl, HybridQuantity::Scalar> const& operand,
    Field<NdArrayImpl, HybridQuantity::Scalar>& derivative) const
{
    this->deriv1D_(operand, derivative);
}




template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::momentsToEx() const
{
    return this->momentsToEx_;
}


template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::momentsToEy() const
{
    return this->momentsToEy_;
}


template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::momentsToEz() const
{
    return this->momentsToEz_;
}

template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::ByToEx() const
{
    return this->ByToEx_;
}

template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::ByToEz() const
{
    return this->ByToEz_;
}

template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BxToEy() const
{
    return this->BxToEy_;
}

template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BxToEz() const
{
    return this->BxToEz_;
}

template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BzToEy() const
{
    return this->BzToEy_;
}

template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BzToEx() const
{
    return this->BzToEx_;
}



template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::ExToMoment() const
{
    return this->ExToMoment_;
}


template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::EyToMoment() const
{
    return this->EyToMoment_;
}


template<std::size_t dim>
LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::EzToMoment() const
{
    return this->EzToMoment_;
}
} // namespace PHARE


#endif // PHARE_CORE_GRID_GRIDLAYOUTYEE_H
