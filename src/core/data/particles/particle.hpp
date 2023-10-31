#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_HPP

#include <array>
#include <random>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <iostream>

#include "core/utilities/point/point.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"


namespace PHARE::core
{
template<typename T = float>
struct ParticleDeltaDistribution
{
    template<typename Generator>
    [[nodiscard]] T operator()(Generator& generator)
    {
        return dist(generator);
    }
    std::uniform_real_distribution<T> dist{0, 1. - std::numeric_limits<T>::epsilon()};
};


template<typename Particle>
[[nodiscard]] auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}


template<size_t dim>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static const size_t dimension = dim;

    Particle(double a_weight, double a_charge, std::array<int, dim> cell,
             std::array<double, dim> a_delta, std::array<double, 3> a_v)
        : weight{a_weight}
        , charge{a_charge}
        , iCell{cell}
        , delta{a_delta}
        , v{a_v}
    {
    }

    Particle() = default;

    double weight;
    double charge;

    std::array<int, dim> iCell    = ConstArray<int, dim>();
    std::array<double, dim> delta = ConstArray<double, dim>();
    std::array<double, 3> v       = ConstArray<double, 3>();

    double Ex = 0, Ey = 0, Ez = 0;
    double Bx = 0, By = 0, Bz = 0;

    [[nodiscard]] bool operator==(Particle<dim> const& that) const
    {
        return (this->weight == that.weight) && //
               (this->charge == that.charge) && //
               (this->iCell == that.iCell) &&   //
               (this->delta == that.delta) &&   //
               (this->v == that.v) &&           //
               (this->Ex == that.Ex) &&         //
               (this->Ey == that.Ey) &&         //
               (this->Ez == that.Ez) &&         //
               (this->Bx == that.Bx) &&         //
               (this->By == that.By) &&         //
               (this->Bz == that.Bz);
    }

    template<std::size_t dimension>
    friend std::ostream& operator<<(std::ostream& out, const Particle<dimension>& particle);
};

template<std::size_t dim>
std::ostream& operator<<(std::ostream& out, Particle<dim> const& particle)
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
    out << ", Exyz : " << particle.Ex << "," << particle.Ey << "," << particle.Ez;
    out << ", Bxyz : " << particle.Bx << "," << particle.By << "," << particle.Bz;
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
    std::array<double, dim>& delta;
    std::array<double, 3>& v;
};




template<std::size_t dim, typename T>
inline constexpr auto is_phare_particle_type
    = std::is_same_v<Particle<dim>, T> or std::is_same_v<ParticleView<dim>, T>;


template<std::size_t dim, template<std::size_t> typename ParticleA,
         template<std::size_t> typename ParticleB>
[[nodiscard]] typename std::enable_if_t<is_phare_particle_type<dim, ParticleA<dim>>
                                            and is_phare_particle_type<dim, ParticleB<dim>>,
                                        bool>
operator==(ParticleA<dim> const& particleA, ParticleB<dim> const& particleB)
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

template<size_t dim, template<std::size_t> typename Particle_t>
[[nodiscard]] typename std::enable_if_t<PHARE::core::is_phare_particle_type<dim, Particle_t<dim>>,
                                        PHARE::core::Particle<dim>>
copy(Particle_t<dim> const& from)
{
    return {from.weight, from.charge, from.iCell, from.delta, from.v};
}


} // namespace std


#endif
