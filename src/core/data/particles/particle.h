#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H

#include "core/utilities/point/point.h"

#include <array>
#include <random>
#include <type_traits>




namespace PHARE::core
{
template<typename T = float>
struct ParticleDeltaDistribution
{
    template<typename Generator>
    float operator()(Generator& generator)
    {
        return dist(generator);
    }
    std::uniform_real_distribution<T> dist{0, 1. - std::numeric_limits<T>::epsilon()};
};


template<std::size_t dim, typename Float = double>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");

    using float_type                   = Float;
    static const std::size_t dimension = dim;

    Float weight;
    Float charge;

    std::array<int32_t, dim> iCell = ZeroArray<int32_t, dim>();
    std::array<float, dim> delta   = ZeroArray<float, dim>();
    std::array<Float, 3> v         = ZeroArray<Float, 3>();

    Float Ex = 0, Ey = 0, Ez = 0;
    Float Bx = 0, By = 0, Bz = 0;
};



template<typename Particle>
auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}



template<std::size_t dim, typename Float = double>
struct ContiguousParticles
{
    std::vector<int> iCell;
    std::vector<float> delta;
    std::vector<Float> weight, charge, v;
    size_t size_;
    auto size() const { return size_; }
    ContiguousParticles(size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * 3)
        , size_(s)
    {
    }
};


} // namespace PHARE::core

#endif
