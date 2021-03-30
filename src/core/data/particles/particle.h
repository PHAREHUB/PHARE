#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H

#include <array>
#include <random>
#include <type_traits>

#include "core/utilities/point/point.h"
#include "core/utilities/span.h"
#include "core/utilities/types.h"


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


template<typename Particle>
auto cellAsPoint(Particle const& particle)
{
    return Point<int, Particle::dimension>{particle.iCell};
}



template<size_t dim>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static const size_t dimension = dim;

    double weight;
    double charge;

    std::array<int, dim> iCell    = ConstArray<int, dim>();
    std::array<double, dim> delta = ConstArray<double, dim>();
    std::array<double, 3> v       = ConstArray<double, 3>();

    double Ex = 0, Ey = 0, Ez = 0;
    double Bx = 0, By = 0, Bz = 0;
};


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



template<std::size_t dim, bool OwnedState = true>
struct ContiguousParticles
{
    static constexpr bool is_contiguous    = true;
    static constexpr std::size_t dimension = dim;
    using ContiguousParticles_             = ContiguousParticles<dim, OwnedState>;

    template<typename T>
    using container_t = std::conditional_t<OwnedState, std::vector<T>, Span<T>>;

    template<bool OS = OwnedState, typename = std::enable_if_t<OS>>
    ContiguousParticles(std::size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * 3)
    {
    }

    template<typename Container_int, typename Container_double>
    ContiguousParticles(Container_int&& _iCell, Container_double&& _delta,
                        Container_double&& _weight, Container_double&& _charge,
                        Container_double&& _v)
        : iCell{_iCell}
        , delta{_delta}
        , weight{_weight}
        , charge{_charge}
        , v{_v}
    {
    }

    std::size_t size() const { return weight.size(); }

    template<std::size_t S, typename T>
    static std::array<T, S>* _array_cast(T const* array)
    {
        return reinterpret_cast<std::array<T, S>*>(const_cast<T*>(array));
    }

    template<typename Return>
    Return _to(std::size_t i)
    {
        return {
            *const_cast<double*>(weight.data() + i),     //
            *const_cast<double*>(charge.data() + i),     //
            *_array_cast<dim>(iCell.data() + (dim * i)), //
            *_array_cast<dim>(delta.data() + (dim * i)), //
            *_array_cast<3>(v.data() + (3 * i)),
        };
    }

    auto copy(std::size_t i) { return _to<Particle<dim>>(i); }
    auto view(std::size_t i) { return _to<ParticleView<dim>>(i); }

    auto operator[](std::size_t i) const { return view(i); }
    auto operator[](std::size_t i) { return view(i); }

    struct iterator
    {
        iterator(ContiguousParticles_* particles)
        {
            for (std::size_t i = 0; i < particles->size(); i++)
                views.emplace_back((*particles)[i]);
        }

        iterator& operator++()
        {
            ++curr_pos;
            return *this;
        }

        bool operator!=(iterator const& other) const { return curr_pos != views.size(); }
        auto& operator*() { return views[curr_pos]; }
        auto& operator*() const { return views[curr_pos]; }

        std::size_t curr_pos = 0;
        std::vector<ParticleView<dim>> views;
    };

    auto begin() { return iterator(this); }
    auto cbegin() const { return iterator(this); }

    auto end() { return iterator(this); }
    auto cend() const { return iterator(this); }

    container_t<int> iCell;
    container_t<double> delta;
    container_t<double> weight, charge, v;
};


template<std::size_t dim>
using ContiguousParticlesView = ContiguousParticles<dim, /*OwnedState=*/false>;



template<std::size_t dim, typename T>
inline constexpr auto is_phare_particle_type
    = std::is_same_v<Particle<dim>, T> or std::is_same_v<ParticleView<dim>, T>;


template<std::size_t dim, template<std::size_t> typename ParticleA,
         template<std::size_t> typename ParticleB>
typename std::enable_if_t<
    is_phare_particle_type<dim, ParticleA<dim>> and is_phare_particle_type<dim, ParticleB<dim>>,
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
typename std::enable_if_t<PHARE::core::is_phare_particle_type<dim, Particle_t<dim>>,
                          PHARE::core::Particle<dim>>
copy(Particle_t<dim> const& from)
{
    return {from.weight, from.charge, from.iCell, from.delta, from.v};
}

} // namespace std


#endif
