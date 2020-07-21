#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H

#include "core/utilities/point/point.h"

#include <array>
#include <random>
#include <type_traits>

#include "kul/tuple.hpp"


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

template<std::size_t dim>
struct Particle
{
};



template<>
struct Particle<1>
{
    double weight;
    double charge;

    std::array<int, 1> iCell   = {{0}};
    std::array<float, 1> delta = {{0.0f}};
    std::array<double, 3> v    = {{0., 0., 0.}};

    double Ex = 0, Ey = 0, Ez = 0;
    double Bx = 0, By = 0, Bz = 0;

    static constexpr std::size_t dimension = 1;
};




template<>
struct Particle<2>
{
    double weight;
    double charge;

    std::array<int, 2> iCell   = {{0, 0}};
    std::array<float, 2> delta = {{0.0f, 0.0f}};
    std::array<double, 3> v    = {{0., 0., 0.}};

    double Ex = 0, Ey = 0, Ez = 0;
    double Bx = 0, By = 0, Bz = 0;

    static constexpr std::size_t dimension = 2;
};



template<>
struct Particle<3>
{
    double weight;
    double charge;

    std::array<int, 3> iCell   = {{0, 0, 0}};
    std::array<float, 3> delta = {{0.f, 0.f, 0.f}};
    std::array<double, 3> v    = {{0., 0., 0.}};

    double Ex = 0, Ey = 0, Ez = 0;
    double Bx = 0, By = 0, Bz = 0;

    static constexpr std::size_t dimension = 3;
};




template<typename Particle>
auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}

template<size_t dim>
struct ContiguousParticleView
{
    static constexpr size_t THREE     = 3;
    static constexpr size_t dimension = dim;

    double& weight;
    double& charge;
    std::array<int, dim>& iCell;
    std::array<float, dim>& delta;
    std::array<double, THREE>& v;
};


template<std::size_t dim, bool OwnedState = true>
struct ContiguousParticles
{
    static constexpr size_t THREE     = 3;
    static constexpr size_t dimension = dim;

    template<typename T>
    using container_t = std::conditional_t<OwnedState, std::vector<T>, kul::Pointers<T>>;

    template<bool OS = OwnedState, typename = typename std::enable_if_t<OS>>
    ContiguousParticles(size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * 3)
    {
    }


    template<typename Container_int, typename Container_float, typename Container_double>
    ContiguousParticles(Container_int&& _iCell, Container_float&& _delta,
                        Container_double&& _weight, Container_double&& _charge,
                        Container_double&& _v)
        : iCell{_iCell}
        , delta{_delta}
        , weight{_weight}
        , charge{_charge}
        , v{_v}
    {
    }

    auto size() const { return weight.size(); }

    template<size_t S, typename T>
    static auto _array_cast(T const* array)
    {
        return reinterpret_cast<std::array<T, S>*>(const_cast<T*>(array));
    }

    auto copy_to_particle(size_t i)
    {
        return Particle<dim>{
            *(weight.data() + i),                        //
            *(charge.data() + i),                        //
            *_array_cast<dim>(iCell.data() + (dim * i)), //
            *_array_cast<dim>(delta.data() + (dim * i)), //
            *_array_cast<THREE>(v.data() + (THREE * i)),
        };
    }

    auto particle(size_t i)
    {
        return ContiguousParticleView<dim>{
            *const_cast<double*>(weight.data() + i),     //
            *const_cast<double*>(charge.data() + i),     //
            *_array_cast<dim>(iCell.data() + (dim * i)), //
            *_array_cast<dim>(delta.data() + (dim * i)), //
            *_array_cast<THREE>(v.data() + (THREE * i)),
        };
    }

    auto operator[](size_t i) { return particle(i); }

    container_t<int> iCell;
    container_t<float> delta;
    container_t<double> weight, charge, v;
};


} // namespace PHARE::core

#endif
