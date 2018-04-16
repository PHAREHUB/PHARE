#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_H
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_H

#include <array>
#include <cstddef>

#include <data/electromag/electromag.h>
#include <data/field/field.h>
#include <data/grid/gridlayout.h>
#include <data/particles/particle.h>
#include <data/vecfield/vecfield.h>

namespace PHARE
{
//! return the size of the index and weights arrays for
//! interpolation at a given order
constexpr int nbrPointsSupport(int order)
{
    return order + 1;
}


//! return the number of points on which to interpolate at order InterpOrder
template<std::size_t InterpOrder>
int computeStartIndex(double normalizedPos)
{
    return static_cast<int>(normalizedPos - (static_cast<double>(InterpOrder) - 1.) / 2.);
}




/** \brief the class Weight aims at computing the weight coefficient for
 *  interpolation at a specific order
 *
 *  This class assumes the interpolation order is known at compile-time
 *  thus there are three specialization for orders 1, 2 and 3.
 *
 *  the class has only one method called computeWeight that takes three arguments:
 *
 *  \param[in] normalized position  is the particle position normalized by grid spacing
 *  \param[in] startIndex first grid index where to interpolate the field
 *  \param[out] weights contains the order+1 weights calculated
 *
 *  these three parameters are given for a specific direction (x, y or z)
 */
template<std::size_t InterpOrder>
class Weighter
{
};



/** \brief Specialization of Weight for first order interpolation
 */
template<>
class Weighter<1>
{
public:
    inline void computeWeight(double normalizedPos, int startIndex,
                              std::array<double, nbrPointsSupport(1)>& weights)
    {
        weights[1] = normalizedPos - static_cast<double>(startIndex);
        weights[0] = 1. - weights[1];
    }

    static const int interp_order = 1;
};


/** \brief specialization of Weighter for second order interpolation
 */
template<>
class Weighter<2>
{
public:
    inline void computeWeight(double normalizedPos, int startIndex,
                              std::array<double, nbrPointsSupport(2)>& weights)
    {
        auto index = startIndex + 1;
        auto delta = static_cast<double>(index) - normalizedPos;
        double coef1, coef2, coef3;
        coef1 = 0.5 + delta;
        coef2 = delta;
        coef3 = 0.5 - delta;

        weights[0] = 0.5 * coef1 * coef1;
        weights[1] = 0.75 - coef2 * coef2;
        weights[2] = 0.5 * coef3 * coef3;
    }

    static const int interp_order = 2;
};



/** \brief specialization of Weighter for third order interpolation
 */
template<>
class Weighter<3>
{
public:
    inline void computeWeight(double normalizedPos, int startIndex,
                              std::array<double, nbrPointsSupport(3)>& weights)
    {
        double coef1, coef2, coef3, coef4;
        auto index = static_cast<double>(startIndex) - normalizedPos;
        coef1      = 1. + 0.5 * index;
        coef2      = index + 1;
        coef3      = index + 2;
        coef4      = 1. - 0.5 * (index + 3);

        double coef2_sq  = coef2 * coef2;
        double coef2_cub = coef2_sq * coef2;
        double coef3_sq  = coef3 * coef3;
        double coef3_cub = coef3_sq * coef3;

        weights[0] = (4. / 3.) * coef1 * coef1 * coef1;
        weights[1] = 2. / 3. - coef2_sq - 0.5 * coef2_cub;
        weights[2] = 2. / 3. - coef3_sq + 0.5 * coef3_cub;
        weights[3] = (4. / 3.) * coef4 * coef4 * coef4;
    }

    static const int interp_order = 3;
};




//! Interpol performs the interpolation of a field using precomputed weights at
//! indices starting at startIndex. The class is templated by the Dimensionality
template<std::size_t dim>
class Interpol
{
};



/** \brief specialization of Interpol for 1D interpolation
 */
template<>
class Interpol<1>
{
public:
    /** Performs the 1D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field
     * \param[in] startIndex is the first of the order+1 indices where to interpolate the field
     * \param[in] weights are the order+1 weights used for the interpolation
     */
    template<typename NdArrayImpl, typename PhysicalQuantity, typename Array1, typename Array2>
    inline double operator()(Field<NdArrayImpl, PhysicalQuantity> const& field,
                             std::array<Centering, 1> const& fieldCentering,
                             Array1 const& startIndex, Array2 const& weights)
    {
        auto particleField      = 0.;
        auto const& xStartIndex = startIndex[static_cast<int>(fieldCentering[0])][0];
        auto const& xWeights    = weights[static_cast<int>(fieldCentering[0])][0];
        auto order_size         = xWeights.size();

        for (auto ik = 0u; ik < order_size; ++ik)
        {
            particleField += field(xStartIndex + ik) * xWeights[ik];
        }
        return particleField;
    }
};


/**\brief Specialization of Interpol for 2D interpolation
 */
template<>
class Interpol<2>
{
public:
    /** Performs the 2D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field in each direction
     * \param[in] startIndex is the first of the order+1 indices where to interpolate the field in
     * both directions \param[in] weights are the order+1 weights used for the interpolation in both
     * directions
     */
    template<typename NdArrayImpl, typename PhysicalQuantity, typename Array1, typename Array2>
    inline double operator()(Field<NdArrayImpl, PhysicalQuantity> const& field,
                             std::array<Centering, 2> const& fieldCentering,
                             Array1 const& startIndex, Array2 const& weights)
    {
        auto const& xStartIndex = startIndex[static_cast<int>(fieldCentering[0])][0];
        auto const& yStartIndex = startIndex[static_cast<int>(fieldCentering[1])][1];
        auto const& xWeights    = weights[static_cast<int>(fieldCentering[0])][0];
        auto const& yWeights    = weights[static_cast<int>(fieldCentering[1])][1];

        auto order_size        = xWeights.size();
        double fieldAtParticle = 0.;
        for (auto ix = 0u; ix < order_size; ++ix)
        {
            double Yinterp = 0.;
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                Yinterp += field(xStartIndex + ix, yStartIndex + iy) * yWeights[iy];
            }
            fieldAtParticle += Yinterp * xWeights[ix];
        }

        return fieldAtParticle;
    }
};



/** \brief Specialization of Interpol for 3D interpolation
 */
template<>
class Interpol<3>
{
public:
    /** Performs the 3D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field in each direction
     * \param[in] startIndex is the first of the order+1 indices where to interpolate the field in
     * the 3 directions \param[in] weights are the order+1 weights used for the interpolation in the
     * 3 directions
     */
    template<typename NdArrayImpl, typename PhysicalQuantity, typename Array1, typename Array2>
    inline double operator()(Field<NdArrayImpl, PhysicalQuantity> const& field,
                             std::array<Centering, 3> const& fieldCentering,
                             Array1 const& startIndex, Array2 const& weights)
    {
        auto const& xStartIndex = startIndex[static_cast<std::size_t>(fieldCentering[0])][0];
        auto const& yStartIndex = startIndex[static_cast<std::size_t>(fieldCentering[1])][1];
        auto const& zStartIndex = startIndex[static_cast<std::size_t>(fieldCentering[2])][2];
        auto const& xWeights    = weights[static_cast<std::size_t>(fieldCentering[0])][0];
        auto const& yWeights    = weights[static_cast<std::size_t>(fieldCentering[1])][1];
        auto const& zWeights    = weights[static_cast<std::size_t>(fieldCentering[2])][2];

        auto order_size        = xWeights.size();
        double fieldAtParticle = 0.;
        for (auto ix = 0u; ix < order_size; ++ix)
        {
            double Yinterp = 0.;
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                double Zinterp = 0.;
                for (auto iz = 0u; iz < order_size; ++iz)
                {
                    Zinterp += field(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)
                               * zWeights[iz];
                }
                Yinterp += Zinterp * yWeights[iy];
            }
            fieldAtParticle += Yinterp * xWeights[ix];
        }
        return fieldAtParticle;
    }
};



/** \brief Interpolator is used to perform particle-mesh interpolations using
 * 1st, 2nd or 3rd order interpolation in 1D, 2D or 3D, on a given layout.
 */
template<std::size_t dim, std::size_t InterpOrder, LayoutType LT>
class Interpolator : private Weighter<InterpOrder>
{
public:
    /**\brief interpolate electromagnetic fields on all particles in the range
     *
     * For each particle :
     *  - The function first calculates the startIndex and weights for interpolation at
     * order InterpOrder and in dimension dim for dual and primal nodes
     *  - then it uses Interpol<> to calculate the interpolation of E and B components
     * onto the particle.
     */
    template<typename PartIterator, typename Electromag>
    inline void operator()(PartIterator begin, PartIterator end, Electromag const& Em)
    {
        // this lambda calculates the startIndex and the order+1 weights for
        // dual field interpolation and puts this at the corresponding location
        // in 'startIndex' and 'weights'. For dual fields, the normalizedPosition
        // is offseted compared to primal ones.
        auto indexAndWeightDual = [this](Particle<dim> const& part) {
            for (auto iDim = 0u; iDim < dim; ++iDim)
            {
                double normalizedPos
                    = part.iCell[iDim] + part.delta[iDim] + dualOffset(InterpOrder);

                startIndex[centering2int(Centering::dual)][iDim]
                    = computeStartIndex<InterpOrder>(normalizedPos);

                weightComputer_.computeWeight(normalizedPos,
                                              startIndex[centering2int(Centering::dual)][iDim],
                                              weights[centering2int(Centering::dual)][iDim]);
            }
        };

        // does the same as above but for a primal field
        auto indexAndWeightPrimal = [this](Particle<dim> const& part) {
            for (auto iDim = 0u; iDim < dim; ++iDim)
            {
                double normalizedPos = part.iCell[iDim] + part.delta[iDim];

                startIndex[centering2int(Centering::primal)][iDim]
                    = computeStartIndex<InterpOrder>(normalizedPos);

                weightComputer_.computeWeight(normalizedPos,
                                              startIndex[centering2int(Centering::primal)][iDim],
                                              weights[centering2int(Centering::primal)][iDim]);
            }
        };

        auto const& Ex = Em.E.getComponent(Component::X);
        auto const& Ey = Em.E.getComponent(Component::Y);
        auto const& Ez = Em.E.getComponent(Component::Z);
        auto const& Bx = Em.B.getComponent(Component::X);
        auto const& By = Em.B.getComponent(Component::Y);
        auto const& Bz = Em.B.getComponent(Component::Z);

        auto const ExCentering = GridLayout<LT, dim>::centering(HybridQuantity::Quantity::Ex);
        auto const EyCentering = GridLayout<LT, dim>::centering(HybridQuantity::Quantity::Ey);
        auto const EzCentering = GridLayout<LT, dim>::centering(HybridQuantity::Quantity::Ez);
        auto const BxCentering = GridLayout<LT, dim>::centering(HybridQuantity::Quantity::Bx);
        auto const ByCentering = GridLayout<LT, dim>::centering(HybridQuantity::Quantity::By);
        auto const BzCentering = GridLayout<LT, dim>::centering(HybridQuantity::Quantity::Bz);


        // for each particle, first calculate the startIndex and weights
        // for dual and primal quantities.
        // then, knowing the centering (primal or dual) of each electromagnetic
        // component, we use Interpol to actually perform the interpolation.
        // the trick here is that the StartIndex and weights have only been calculated
        // twice, and not for each E,B component.
        for (auto currPart = begin; currPart != end; ++currPart)
        {
            indexAndWeightPrimal(*currPart);
            indexAndWeightDual(*currPart);

            currPart->Ex = interpol_(Ex, ExCentering, startIndex, weights);
            currPart->Ey = interpol_(Ey, EyCentering, startIndex, weights);
            currPart->Ez = interpol_(Ez, EzCentering, startIndex, weights);
            currPart->Bx = interpol_(Bx, BxCentering, startIndex, weights);
            currPart->By = interpol_(By, ByCentering, startIndex, weights);
            currPart->Bz = interpol_(Bz, BzCentering, startIndex, weights);
        }
    }


private:
    static_assert(dim <= 3 && dim > 0 && InterpOrder >= 1 && InterpOrder <= 3, "error");

    Weighter<InterpOrder> weightComputer_;
    Interpol<dim> interpol_;

    // array[dual/primal][dim]
    std::array<std::array<int, dim>, 2> startIndex;
    std::array<std::array<std::array<double, nbrPointsSupport(InterpOrder)>, dim>, 2> weights;

    /**
     * @brief dualOffset returns the offset by which changing the
     * startIndex for dual node interpolation. This offset depends on
     * interpolation order. */
    inline constexpr double dualOffset(int order)
    {
        std::array<double, 3> offsets = {{-0.5, -0.5, 0.5}};
        return offsets[static_cast<std::array<double, 3>::size_type>(order - 1)];
    }
};




} // namespace PHARE

#endif
