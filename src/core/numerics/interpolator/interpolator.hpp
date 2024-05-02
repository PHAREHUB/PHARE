#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_HPP



#include <array>
#include <cstddef>

#include "core/utilities/point/point.hpp"
#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/range/range.hpp"



namespace PHARE::core
{


//! return the size of the index and weights arrays for
//! interpolation at a given order
// Number of points where the interpolator of order interpOrder
// will deposit mass/momentum. This is hence the size of the
// index and weight arrays for interpolation at a given order
constexpr int nbrPointsSupport(int interpOrder)
{
    return interpOrder + 1;
}



/** \brief the class Weight aims at computing the weight coefficient for
 *  interpolation at a specific order
 *
 *  This class assumes the interpolation order is known at compile-time
 *  thus there are three specialization for orders 1, 2 and 3.
 *
 *  the class has only one method called computeWeight that takes three arguments:
 *
 *  \param[in] normalized position is the particle position in the cell normalized by grid
 * spacing \param[in] startIndex first grid index where to interpolate the field \param[out]
 * weights contains the nbrPointsSupport weights calculated
 *
 *  these three parameters are given for a specific direction (x, y or z)
 */
template<std::size_t interpOrder>
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

    static constexpr int interp_order = 1;
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

    static constexpr int interp_order = 2;
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
        constexpr double _4_over_3 = 4. / 3.;
        constexpr double _2_over_3 = 2. / 3.;

        auto index   = static_cast<double>(startIndex) - normalizedPos;
        double coef1 = 1. + 0.5 * index;
        double coef2 = index + 1;
        double coef3 = index + 2;
        double coef4 = 1. - 0.5 * (index + 3);

        double coef2_sq  = coef2 * coef2;
        double coef2_cub = coef2_sq * coef2;
        double coef3_sq  = coef3 * coef3;
        double coef3_cub = coef3_sq * coef3;

        weights[0] = _4_over_3 * coef1 * coef1 * coef1;
        weights[1] = _2_over_3 - coef2_sq - 0.5 * coef2_cub;
        weights[2] = _2_over_3 - coef3_sq + 0.5 * coef3_cub;
        weights[3] = _4_over_3 * coef4 * coef4 * coef4;
    }

    static constexpr int interp_order = 3;
};




//! Interpol performs the interpolation of a field using precomputed weights at
//! indices starting at startIndex. The class is templated by the Dimensionality
template<std::size_t dim>
class MeshToParticle
{
};


template<std::size_t dimdex, typename GridLayout, auto quantity, typename IndexWeights>
NO_DISCARD auto static start_index_and_weights_for_qty(IndexWeights const& indexWeights)
{
    auto constexpr centerings                              = GridLayout::centering(quantity);
    auto const& [d_starts, d_weights, p_starts, p_weights] = indexWeights;

    if constexpr (centerings[dimdex] == QtyCentering::primal)
        return std::forward_as_tuple(p_starts[dimdex], p_weights[dimdex]);
    else
        return std::forward_as_tuple(d_starts[dimdex], d_weights[dimdex]);
}

/** \brief specialization of Interpol for 1D interpolation
 */
template<>
class MeshToParticle<1>
{
public:
    /** Performs the 1D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field
     * \param[in] startIndex is the first of the nbrPointsSupport indices where to interpolate
     * the field \param[in] weights are the nbrPointsSupport weights used for the interpolation
     */
    template<typename GridLayout, auto quantity, typename Field, typename IndexWeights>
    NO_DISCARD inline auto operator()(Field const& field, IndexWeights const& indexWeights)
    {
        auto const& [xStartIndex, xWeights]
            = start_index_and_weights_for_qty<0, GridLayout, quantity>(indexWeights);

        auto const& order_size = xWeights.size();
        auto fieldAtParticle   = 0.;

        for (auto ik = 0u; ik < order_size; ++ik)
        {
            fieldAtParticle += field(xStartIndex + ik) * xWeights[ik];
        }
        return fieldAtParticle;
    }
};


/**\brief Specialization of Interpol for 2D interpolation
 */
template<>
class MeshToParticle<2>
{
public:
    /** Performs the 2D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field in each
     * direction \param[in] startIndex is the first of the nbrPointsSupport indices where to
     * interpolate the field in both directions \param[in] weights are the nbrPointsSupport
     * weights used for the interpolation in both directions
     */
    template<typename GridLayout, auto quantity, typename Field, typename IndexWeights>
    NO_DISCARD inline auto operator()(Field const& field, IndexWeights const& indexWeights)
    {
        auto const& [xStartIndex, xWeights]
            = start_index_and_weights_for_qty<0, GridLayout, quantity>(indexWeights);
        auto const& [yStartIndex, yWeights]
            = start_index_and_weights_for_qty<1, GridLayout, quantity>(indexWeights);

        auto const& order_size = xWeights.size();

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
class MeshToParticle<3>
{
public:
    /** Performs the 3D interpolation
     * \param[in] field is the field from which values are interpolated
     * \param[in] fieldCentering is the centering (dual or primal) of the field in each
     * direction \param[in] startIndex is the first of the nbrPointsSupport indices where to
     * interpolate the field in the 3 directions \param[in] weights are the nbrPointsSupport
     * weights used for the interpolation in the 3 directions
     */
    template<typename GridLayout, auto quantity, typename Field, typename IndexWeights>
    NO_DISCARD inline auto operator()(Field const& field, IndexWeights const& indexWeights)
    {
        auto const& [xStartIndex, xWeights]
            = start_index_and_weights_for_qty<0, GridLayout, quantity>(indexWeights);
        auto const& [yStartIndex, yWeights]
            = start_index_and_weights_for_qty<1, GridLayout, quantity>(indexWeights);
        auto const& [zStartIndex, zWeights]
            = start_index_and_weights_for_qty<2, GridLayout, quantity>(indexWeights);

        auto const& order_size = xWeights.size();

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




//! ParticleToMesh projects a particle density and flux to given grids
template<std::size_t dim>
class ParticleToMesh
{
};



/** \brief specialization of ParticleToMesh for 1D interpolation */
template<>
class ParticleToMesh<1>
{
public:
    template<typename Field, typename Particle, typename Func, typename Indexes, typename Weights>
    inline void operator()(Field& field, Particle const& particle, Func&& func,
                           Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xStartIndex] = startIndex;
        auto const& [xWeights]    = weights;
        auto const& order_size    = xWeights.size();

        auto const deposit = func(particle) * particle.weight * coef;

        for (auto ik = 0u; ik < order_size; ++ik)
        {
            field(xStartIndex + ik) += deposit * xWeights[ik];
        }
    }
};




/** \brief specialization of ParticleToMesh for 2D interpolation */
template<>
class ParticleToMesh<2>
{
public:
    template<typename Field, typename Particle, typename Func, typename Indexes, typename Weights>
    inline void operator()(Field& field, Particle const& particle, Func&& func,
                           Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xStartIndex, yStartIndex] = startIndex;
        auto const& [xWeights, yWeights]       = weights;
        auto const& order_size                 = xWeights.size();

        auto const deposit = func(particle) * particle.weight * coef;

        for (auto ix = 0u; ix < order_size; ++ix)
        {
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                auto x = xStartIndex + ix;
                auto y = yStartIndex + iy;

                field(x, y) += deposit * xWeights[ix] * yWeights[iy];
            }
        }
    }
};




/** \brief specialization of ParticleToMesh for 3D interpolation */
template<>
class ParticleToMesh<3>
{
public:
    template<typename Field, typename Particle, typename Func, typename Indexes, typename Weights>
    inline void operator()(Field& field, Particle const& particle, Func&& func,
                           Indexes const& startIndex, Weights const& weights, double coef = 1.)
    {
        auto const& [xStartIndex, yStartIndex, zStartIndex] = startIndex;
        auto const& [xWeights, yWeights, zWeights]          = weights;
        auto const& order_size                              = xWeights.size();

        auto const deposit = func(particle) * particle.weight * coef;

        for (auto ix = 0u; ix < order_size; ++ix)
        {
            for (auto iy = 0u; iy < order_size; ++iy)
            {
                for (auto iz = 0u; iz < order_size; ++iz)
                {
                    auto x = xStartIndex + ix;
                    auto y = yStartIndex + iy;
                    auto z = zStartIndex + iz;

                    field(x, y, z) += deposit * xWeights[ix] * yWeights[iy] * zWeights[iz];
                }
            }
        }
    }
};




/** \brief Interpolator is used to perform particle-mesh interpolations using
 * 1st, 2nd or 3rd order interpolation in 1D, 2D or 3D, on a given layout.
 */
template<std::size_t dim, std::size_t interpOrder>
class Interpolator : private Weighter<interpOrder>
{
protected:
    // this calculates the startIndex and the nbrPointsSupport() weights for
    // dual field interpolation and puts this at the corresponding location
    // in 'startIndex' and 'weights'. For dual fields, the normalizedPosition
    // is offseted compared to primal ones.
    template<typename CenteringT, CenteringT centering, typename GridLayout, typename ICell,
             typename Delta>
    auto indexAndWeights_(GridLayout const& layout, ICell const& iCell_, Delta const& delta)
    {
        // dual weights require -.5 to take the correct position weight
        auto constexpr dual_offset = .5;

        auto const& [startIndex_, weights_] = [&]() {
            if constexpr (centering == QtyCentering::dual)
                return std::forward_as_tuple(dual_startIndex_, dual_weights_);
            else
                return std::forward_as_tuple(primal_startIndex_, primal_weights_);
        }();

        auto iCell = layout.AMRToLocal(Point{iCell_});
        for (auto iDim = 0u; iDim < dimension; ++iDim)
        {
            startIndex_[iDim]
                = iCell[iDim] - computeStartLeftShift<CenteringT, centering>(delta[iDim]);

            double normalizedPos = iCell[iDim] + delta[iDim];

            if constexpr (centering == QtyCentering::dual)
                normalizedPos -= dual_offset;

            weightComputer_.computeWeight(normalizedPos, startIndex_[iDim], weights_[iDim]);
        }
    }

public:
    auto static constexpr interp_order = interpOrder;
    auto static constexpr dimension    = dim;

    /**\brief interpolate electromagnetic fields on a particle and return particles EB
     *
     * For each particle :
     *  - The function first calculates the startIndex and weights for interpolation at
     * order InterpOrder and in dimension dim for dual and primal nodes
     *  - then it uses Interpol<> to calculate the interpolation of E and B components
     * onto the particle.
     */
    template<typename Particle_t, typename Electromag, typename GridLayout>
    inline auto operator()(Particle_t& currPart, Electromag const& Em, GridLayout const& layout)
    {
        using E_B_tuple = std::tuple<std::array<double, 3>, std::array<double, 3>>;
        using Scalar    = HybridQuantity::Scalar;

        // for each particle, first calculate the startIndex and weights for dual and
        // primal quantities. then, knowing the centering (primal or dual) of each
        // electromagnetic component, we use Interpol to actually perform the
        // interpolation. the trick here is that the StartIndex and weights have only been
        // calculated twice, and not for each E,B component.

        auto& iCell = currPart.iCell;
        auto& delta = currPart.delta;
        indexAndWeights_<QtyCentering, QtyCentering::dual>(layout, iCell, delta);
        indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, iCell, delta);

        auto indexWeights = std::forward_as_tuple(dual_startIndex_, dual_weights_,
                                                  primal_startIndex_, primal_weights_);

        E_B_tuple particle_EB;
        auto& [pE, pB]        = particle_EB;
        auto& [pEx, pEy, pEz] = pE;
        auto& [pBx, pBy, pBz] = pB;

        auto const& [Ex, Ey, Ez] = Em.E();
        auto const& [Bx, By, Bz] = Em.B();

        pEx = meshToParticle_.template operator()<GridLayout, Scalar::Ex>(Ex, indexWeights);
        pEy = meshToParticle_.template operator()<GridLayout, Scalar::Ey>(Ey, indexWeights);
        pEz = meshToParticle_.template operator()<GridLayout, Scalar::Ez>(Ez, indexWeights);
        pBx = meshToParticle_.template operator()<GridLayout, Scalar::Bx>(Bx, indexWeights);
        pBy = meshToParticle_.template operator()<GridLayout, Scalar::By>(By, indexWeights);
        pBz = meshToParticle_.template operator()<GridLayout, Scalar::Bz>(Bz, indexWeights);

        return particle_EB;
    }



    /**\brief interpolate electromagnetic fields on all particles in the range
     *
     * For each particle :
     *  - The function first calculates the startIndex and weights for interpolation at
     * order InterpOrder and in dimension dim for dual and primal nodes
     *  - then it uses Interpol<> to calculate the interpolation of E and B components
     * onto the particle.
     */
    template<typename ParticleRange, typename VecField, typename GridLayout, typename Field>
    inline void operator()(ParticleRange& particleRange, Field& density, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
    {
        auto begin                        = particleRange.begin();
        auto end                          = particleRange.end();
        auto& startIndex_                 = primal_startIndex_;
        auto& weights_                    = primal_weights_;
        auto const& [xFlux, yFlux, zFlux] = flux();


        PHARE_LOG_START(3, "ParticleToMesh::operator()");

        for (auto currPart = begin; currPart != end; ++currPart)
        {
            indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, currPart->iCell,
                                                                 currPart->delta);

            particleToMesh_(
                density, *currPart, [](auto const& part) { return 1.; }, startIndex_, weights_,
                coef);
            particleToMesh_(
                xFlux, *currPart, [](auto const& part) { return part.v[0]; }, startIndex_, weights_,
                coef);
            particleToMesh_(
                yFlux, *currPart, [](auto const& part) { return part.v[1]; }, startIndex_, weights_,
                coef);
            particleToMesh_(
                zFlux, *currPart, [](auto const& part) { return part.v[2]; }, startIndex_, weights_,
                coef);
        }
        PHARE_LOG_STOP(3, "ParticleToMesh::operator()");
    }
    template<typename ParticleRange, typename VecField, typename GridLayout, typename Field>
    inline void operator()(ParticleRange&& range, Field& density, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
    {
        (*this)(range, density, flux, layout, coef);
    }




    /**
     * @brief Given a delta and an interpolation order, deduce which lower index to start
     * traversing from
     */
    template<typename CenteringT, CenteringT Centering>
    NO_DISCARD static int computeStartLeftShift([[maybe_unused]] double delta)
    {
        static_assert(interpOrder > 0 and interpOrder < 4);

        // If this is no longer true, it should be handled here via if constexpr/etc

        if constexpr (interpOrder == 1)
        {
            if constexpr (Centering == QtyCentering::primal)
                return 0;
            else
                return (delta < .5 ? 1 : 0);
        }

        else if constexpr (interpOrder == 2)
        {
            if constexpr (Centering == QtyCentering::primal)
                return (delta < .5 ? 1 : 0);
            else
                return 1;
        }

        else if constexpr (interpOrder == 3)
        {
            if constexpr (Centering == QtyCentering::primal)
                return 1;
            else
                return (delta < .5 ? 2 : 1);
        }
    }


protected:
    static_assert(dimension <= 3 && dimension > 0 && interpOrder >= 1 && interpOrder <= 3, "error");

    using Starts  = std::array<std::uint32_t, dimension>;
    using Weights = std::array<std::array<double, nbrPointsSupport(interpOrder)>, dimension>;

    Weighter<interpOrder> weightComputer_;
    MeshToParticle<dimension> meshToParticle_;
    ParticleToMesh<dimension> particleToMesh_;

    Starts dual_startIndex_;
    Weights dual_weights_;

    Starts primal_startIndex_;
    Weights primal_weights_;
};

template<std::size_t dim, std::size_t interpOrder>
class MomentumTensorInterpolator : public Interpolator<dim, interpOrder>
{
public:
    template<typename ParticleRange, typename TensorField, typename GridLayout>
    inline void operator()(ParticleRange&& particleRange, TensorField& momentumTensor,
                           GridLayout const& layout, double mass = 1.)
    {
        auto begin                           = particleRange.begin();
        auto end                             = particleRange.end();
        auto& startIndex_                    = this->primal_startIndex_;
        auto& weights_                       = this->primal_weights_;
        auto const& [xx, xy, xz, yy, yz, zz] = momentumTensor();

        PHARE_LOG_START(3, "ParticleToMesh::operator()");

        for (auto currPart = begin; currPart != end; ++currPart)
        {
            this->template indexAndWeights_<QtyCentering, QtyCentering::primal>(
                layout, currPart->iCell, currPart->delta);

            this->particleToMesh_(
                xx, *currPart, [](auto const& part) { return part.v[0] * part.v[0]; }, startIndex_,
                weights_, mass);
            this->particleToMesh_(
                xy, *currPart, [](auto const& part) { return part.v[0] * part.v[1]; }, startIndex_,
                weights_, mass);
            this->particleToMesh_(
                xz, *currPart, [](auto const& part) { return part.v[0] * part.v[2]; }, startIndex_,
                weights_, mass);
            this->particleToMesh_(
                yy, *currPart, [](auto const& part) { return part.v[1] * part.v[1]; }, startIndex_,
                weights_, mass);
            this->particleToMesh_(
                yz, *currPart, [](auto const& part) { return part.v[1] * part.v[2]; }, startIndex_,
                weights_, mass);
            this->particleToMesh_(
                zz, *currPart, [](auto const& part) { return part.v[2] * part.v[2]; }, startIndex_,
                weights_, mass);
        }
        PHARE_LOG_STOP(3, "ParticleToMesh::operator()");
    }
};

} // namespace PHARE::core

#endif
