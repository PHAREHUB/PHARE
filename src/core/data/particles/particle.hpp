#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"

#include <array>
#include <random>
#include <iostream>
#include <algorithm>
#include <type_traits>


namespace PHARE::core
{


template<size_t dim>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static auto constexpr dimension = dim;

    constexpr Particle() _PHARE_ALL_FN_ {}
    Particle(Particle const&)            = default;
    Particle(Particle&&)                 = default;
    Particle& operator=(Particle const&) = default;
    Particle& operator=(Particle&&)      = default;

    Particle(double const& a_weight, double const& a_charge, std::array<int, dim> const& cell,
             std::array<double, dim> const& a_delta,
             std::array<double, 3> const& a_v) _PHARE_ALL_FN_ //
        : weight_{a_weight},
          charge_{a_charge},
          iCell_{cell},
          delta_{a_delta},
          v_{a_v}
    {
    }


    template<typename Particle_t>
    Particle(Particle_t const& p) _PHARE_ALL_FN_ //
        : weight_{p.weight()},
          charge_{p.charge()},
          iCell_{p.iCell()},
          delta_{p.delta()},
          v_{p.v()}
    {
    }


    //                                                             1d  2d  3d
    double weight_                 = 0;                         // 8   8   8
    double charge_                 = 0;                         // 16  16  16
    std::array<int, dim> iCell_    = ConstArray<int, dim>();    // 20  24  28
    std::array<double, dim> delta_ = ConstArray<double, dim>(); // 28  40  52
    std::array<double, 3> v_       = ConstArray<double, 3>();   // 52  64  76

    NO_DISCARD bool operator==(Particle<dim> const& that) const
    {
        return (this->weight_ == that.weight_) && //
               (this->charge_ == that.charge_) && //
               (this->iCell_ == that.iCell_) &&   //
               (this->delta_ == that.delta_) &&   //
               (this->v_ == that.v_);
    }

    auto& weight() _PHARE_ALL_FN_ { return weight_; }
    auto& weight() const _PHARE_ALL_FN_ { return weight_; }

    auto& charge() _PHARE_ALL_FN_ { return charge_; }
    auto& charge() const _PHARE_ALL_FN_ { return charge_; }

    auto& iCell() _PHARE_ALL_FN_ { return iCell_; }
    auto& iCell() const _PHARE_ALL_FN_ { return iCell_; }

    auto& delta() _PHARE_ALL_FN_ { return delta_; }
    auto& delta() const _PHARE_ALL_FN_ { return delta_; }

    auto& v() _PHARE_ALL_FN_ { return v_; }
    auto& v() const _PHARE_ALL_FN_ { return v_; }

    template<std::size_t dimension>
    friend std::ostream& operator<<(std::ostream& out, Particle<dimension> const& particle);


    auto copy() const { return *this; }
};




template<typename T = float>
struct ParticleDeltaDistribution
{
    template<typename Generator>
    NO_DISCARD T operator()(Generator& generator)
    { return dist(generator); }
    std::uniform_real_distribution<T> dist{0, 1. - std::numeric_limits<T>::epsilon()};
};


template<typename T, std::size_t dim>
NO_DISCARD auto cellAsPoint(std::array<T, dim> const& iCell) _PHARE_ALL_FN_
{ return Point<int, dim>{iCell}; }


template<typename Particle>
NO_DISCARD auto cellAsPoint(Particle const& particle) _PHARE_ALL_FN_
{ return cellAsPoint(particle.iCell()); }



template<std::size_t dim>
struct ParticlePosition
{
    static size_t const dimension = dim;

    std::array<int, dim> _iCell;
    std::array<double, dim> _delta;

    auto& iCell() const { return _iCell; }
    auto& delta() const { return _delta; }
};



template<size_t dim>
struct CountedParticle : public Particle<dim>
{
    CountedParticle(double const& a_weight, double const& a_charge,
                    std::array<int, dim> const& cell, std::array<double, dim> const& a_delta,
                    std::array<double, 3> const& a_v, std::size_t const& id_ = 0) _PHARE_ALL_FN_
        : Particle<dim>{a_weight, a_charge, cell, a_delta, a_v},
          id{id_}
    {
    }
    CountedParticle() _PHARE_ALL_FN_ {}

    std::size_t id = 0;

    auto copy() const { return *this; }
    auto& super() const { return *this; }

    template<std::size_t dimension>
    friend std::ostream& operator<<(std::ostream& out, CountedParticle<dimension> const& particle);
};

template<size_t dim>
auto to_string(Particle<dim> const& particle)
{
    std::stringstream out;
    auto const append = [&](std::string const name, auto const& value) {
        out << name << "(" << value[0];
        for (std::size_t i = 1; i < dim; ++i)
            out << "," << value[i];
        out << "), ";
    };
    append("iCell", particle.iCell());
    append("delta", particle.delta());
    append("v", particle.v());
    out << "), charge: " << particle.charge() << ", weight: " << particle.weight();
    return out.str();
}

template<template<std::size_t> typename Particle_t, std::size_t dim, typename Stream>
Stream& write_to_stream(Particle_t<dim> const& particle, Stream& out, bool const new_line = true)
{
    out << to_string(particle);
    if constexpr (std::is_same_v<Particle_t<dim>, CountedParticle<dim>>)
        out << ", id : " << particle.id;
    if (new_line)
        out << '\n';
    return out;
}

template<std::size_t dim>
std::ostream& operator<<(std::ostream& out, Particle<dim> const& particle)
{ return write_to_stream(particle, out); }

template<std::size_t dim>
std::ostream& operator<<(std::ostream& out, CountedParticle<dim> const& particle)
{
    write_to_stream(particle, out, /*new_line =*/false);
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


template<typename Box_t, typename RValue = std::size_t>
struct CellFlattener
{
    template<typename Icell>
    NO_DISCARD RValue operator()(Icell const& icell) const
    {
        if constexpr (Box_t::dimension == 2)
            return icell[1] + icell[0] * shape[1];
        if constexpr (Box_t::dimension == 3)
            return icell[2] + icell[1] * shape[2] + icell[0] * shape[1] * shape[2];
        return icell[0];
    }

    Box_t const box;
    std::array<int, Box_t::dimension> shape = box.shape().toArray();
};


template<std::size_t dim, typename T>
inline constexpr auto is_phare_particle_type
    = std::is_same_v<Particle<dim>, T> or std::is_same_v<ParticleView<dim>, T>;


template<std::size_t dim, template<std::size_t> typename ParticleA,
         template<std::size_t> typename ParticleB>
NO_DISCARD typename std::enable_if_t<is_phare_particle_type<dim, ParticleA<dim>>
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

template<std::size_t dim>
void check_particle(Particle<dim> const& p) _PHARE_ALL_FN_
{
#if PHARE_DEBUG && !PHARE_HAVE_GPU
    for (std::size_t i = 0; i < dim; ++i)
    {
        assert(not std::isnan(p.v()[i]));
        assert(not std::isnan(p.delta()[i]));
    }
#endif
}


auto shift_particle(auto particle, auto const& shift)
{
    array_op<PlusEquals<int>>(particle.iCell(), shift);
    return particle;
}


} // namespace PHARE::core


template<std::size_t dim>
void swap(PHARE::core::Particle<dim>& a, PHARE::core::Particle<dim>& b)
{
    // if (a == Particle<dim>{})
    //     std::abort();
    PHARE_LOG_LINE_STR(a);
    // PHARE_LOG_LINE_STR(b);
    if (&a != &b)
        std::swap(a, b);
}

#endif
