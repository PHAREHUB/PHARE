#ifndef PHARE_BENCH_CORE_NUMERICS_INTERPOLATOR_BENCH_INTERPOLATOR_ACCESS_HPP
#define PHARE_BENCH_CORE_NUMERICS_INTERPOLATOR_BENCH_INTERPOLATOR_ACCESS_HPP

#include "core/numerics/interpolator/interpolator.hpp"

#include <array>
#include <cstddef>
#include <tuple>


// frozen copy of the MeshToParticle/ParticleToMesh field access before stride based access
//  and the fused __restrict deposit, as a baseline for the production Interpolator
//  index and weight computation is the production one (Interpolator::indexAndWeights_)

namespace PHARE::core::bench
{

// a field index computed per stencil point

template<std::size_t dim>
class PreviousMeshToParticle
{
};

template<>
class PreviousMeshToParticle<2>
{
public:
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

template<>
class PreviousMeshToParticle<3>
{
public:
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



template<std::size_t dim>
class PreviousParticleToMesh
{
};

template<>
class PreviousParticleToMesh<2>
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

template<>
class PreviousParticleToMesh<3>
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



/** @brief The production Interpolator's index and weight computation, with the previous
 * MeshToParticle/ParticleToMesh: one ParticleToMesh call per field (density and 3 flux).
 */
template<std::size_t dim, std::size_t interpOrder>
class PreviousInterpolator : public Interpolator<dim, interpOrder>
{
    using Super = Interpolator<dim, interpOrder>;
    using Super::dual_startIndex_;
    using Super::dual_weights_;
    using Super::primal_startIndex_;
    using Super::primal_weights_;

public:
    template<typename Particle_t, typename Electromag, typename GridLayout>
    inline auto operator()(Particle_t& currPart, Electromag const& Em, GridLayout const& layout)
    {
        using E_B_tuple = std::tuple<std::array<double, 3>, std::array<double, 3>>;
        using Scalar    = HybridQuantity::Scalar;

        auto& iCell = currPart.iCell;
        auto& delta = currPart.delta;
        this->template indexAndWeights_<QtyCentering, QtyCentering::dual>(layout, iCell, delta);
        this->template indexAndWeights_<QtyCentering, QtyCentering::primal>(layout, iCell, delta);

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

    template<typename ParticleRange, typename VecField, typename GridLayout, typename Field>
    inline void operator()(ParticleRange& particleRange, Field& particleDensity, VecField& flux,
                           GridLayout const& layout, double coef = 1.)
    {
        auto& startIndex_                 = primal_startIndex_;
        auto& weights_                    = primal_weights_;
        auto const& [xFlux, yFlux, zFlux] = flux();

        for (auto currPart = particleRange.begin(); currPart != particleRange.end(); ++currPart)
        {
            this->template indexAndWeights_<QtyCentering, QtyCentering::primal>(
                layout, currPart->iCell, currPart->delta);

            particleToMesh_(
                particleDensity, *currPart, [](auto const& part) { return 1.; }, startIndex_,
                weights_, coef);
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
    }

private:
    PreviousMeshToParticle<dim> meshToParticle_;
    PreviousParticleToMesh<dim> particleToMesh_;
};


// bench versions: 0 = previous (frozen copy), 1 = production Interpolator
template<std::size_t dim, std::size_t interpOrder, int version>
using Interpolator_t = std::tuple_element_t<
    version, std::tuple<PreviousInterpolator<dim, interpOrder>, Interpolator<dim, interpOrder>>>;

} // namespace PHARE::core::bench


#endif /* PHARE_BENCH_CORE_NUMERICS_INTERPOLATOR_BENCH_INTERPOLATOR_ACCESS_HPP */
