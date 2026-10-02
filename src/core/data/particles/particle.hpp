#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_HPP


#include "core/def.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/particles/particle_storage.hpp"


#include <array>
#include <random>
#include <iostream>
#include <algorithm>
#include <type_traits>


namespace PHARE::core
{
template<typename T = float>
struct ParticleDeltaDistribution
{
    template<typename Generator>
    NO_DISCARD T operator()(Generator& generator)
    {
        return dist(generator);
    }
    std::uniform_real_distribution<T> dist{0, 1. - std::numeric_limits<T>::epsilon()};
};


template<typename Particle>
NO_DISCARD auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}


template<size_t dim, typename Delta_t = ParticleDelta_t, typename V_t = ParticleV_t>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static std::size_t constexpr dimension = dim;
    using delta_type                       = Delta_t;
    using v_type                           = V_t;

    Particle(double a_weight, double a_charge, std::array<int, dim> cell,
             std::array<double, dim> a_delta, std::array<double, 3> a_v)
        : weight{a_weight}
        , charge{a_charge}
        , iCell{cell}
        , delta{a_delta}
        , v{a_v}
    {
    }

    template<typename D, typename V>
    Particle(double a_weight, double a_charge, std::array<int, dim> cell,
             MultiPrecisionArray<D, dim> const& a_delta, MultiPrecisionArray<V, 3> const& a_v)
        : weight{a_weight}
        , charge{a_charge}
        , iCell{cell}
        , delta{a_delta}
        , v{a_v}
    {
    }

    Particle() = default;

    double weight = 0;
    double charge = 0;

    // {} zero initialization
    std::array<int, dim> iCell{};
    MultiPrecisionArray<Delta_t, dim> delta{};
    MultiPrecisionArray<V_t, 3> v{};

    NO_DISCARD bool operator==(Particle const& that) const
    {
        return (this->weight == that.weight) && //
               (this->charge == that.charge) && //
               (this->iCell == that.iCell) &&   //
               (this->delta == that.delta) &&   //
               (this->v == that.v);
    }
};

template<std::size_t dim, typename Delta_t, typename V_t>
std::ostream& operator<<(std::ostream& out, Particle<dim, Delta_t, V_t> const& particle)
{
    out << "iCell(";
    for (auto c : particle.iCell)
    {
        out << c << ",";
    }
    out << "), delta(";
    for (auto d : particle.delta)
    {
        out << d << ",";
    }
    out << "), v(";
    for (auto v : particle.v)
    {
        out << v << ",";
    }
    out << "), charge : " << particle.charge << ", weight : " << particle.weight;
    out << '\n';
    return out;
}


template<std::size_t dim>
struct ParticleView
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static constexpr std::size_t dimension = dim;

    double& weight;
    double& charge;
    std::array<int, dim>& iCell;
    MultiPrecisionArray<double, dim>& delta;
    MultiPrecisionArray<double, 3>& v;
};



template<typename T>
struct is_phare_particle : std::false_type
{
};
template<std::size_t dim, typename Delta_t, typename V_t>
struct is_phare_particle<Particle<dim, Delta_t, V_t>> : std::true_type
{
};
template<std::size_t dim>
struct is_phare_particle<ParticleView<dim>> : std::true_type
{
};

template<typename T>
inline constexpr bool is_phare_particle_v = is_phare_particle<T>::value;


template<typename ParticleA, typename ParticleB>
    requires(is_phare_particle_v<ParticleA> and is_phare_particle_v<ParticleB>
             and ParticleA::dimension == ParticleB::dimension)
NO_DISCARD bool operator==(ParticleA const& particleA, ParticleB const& particleB)
{
    return particleA.weight == particleB.weight and //
           particleA.charge == particleB.charge and //
           particleA.iCell == particleB.iCell and   //
           particleA.delta == particleB.delta and   //
           particleA.v == particleB.v;
}

} // namespace PHARE::core


namespace std
{

template<size_t dim, typename Delta_t, typename V_t>
NO_DISCARD PHARE::core::Particle<dim, Delta_t, V_t>
copy(PHARE::core::Particle<dim, Delta_t, V_t> const& from)
{
    return from;
}

template<size_t dim>
NO_DISCARD PHARE::core::Particle<dim> copy(PHARE::core::ParticleView<dim> const& from)
{
    return {from.weight, from.charge, from.iCell, from.delta, from.v};
}


} // namespace std


#endif
