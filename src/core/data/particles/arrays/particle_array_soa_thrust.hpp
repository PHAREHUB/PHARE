#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP


#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"


#if PHARE_HAVE_THRUST
#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>
#include <thrust/iterator/zip_iterator.h>
#endif // PHARE_HAVE_THRUST

#include <type_traits>


namespace PHARE::core
{

// used in SOA_ParticelArray::operator[](std::size_t) for std algorithm interactions
//  allows copying between vector and array impls without iterators
//   could be a better way maybe but it's not so simple to change
template<std::size_t dim>
using SoAParticle_crt = std::tuple<double const&,                  //  weight
                                   double const&,                  // charge
                                   std::array<int, dim> const&,    // iCell
                                   std::array<double, dim> const&, // delta
                                   std::array<double, 3> const&    // v

                                   >;
template<std::size_t dim>
using SoAParticle_rt = std::tuple<double&,                  //  weight
                                  double&,                  // charge
                                  std::array<int, dim>&,    // iCell
                                  std::array<double, dim>&, // delta
                                  std::array<double, 3>&    // v
                                  >;




#if PHARE_HAVE_THRUST
template<typename SoAParticles_t>
struct SoAZipParticle;


struct SoAIteratorAdaptor
{
    template<typename Particles>
    auto static make(Particles& particles, std::size_t const& idx = 0) _PHARE_ALL_FN_
    {
        if constexpr (Particles::storage_mode == StorageMode::SPAN)
        {
            return thrust::make_zip_iterator( //
                thrust::make_tuple(           //
                    particles.weight_ + idx,  // 0
                    particles.charge_ + idx,  // 1
                    particles.iCell_ + idx,   // 2
                    particles.delta_ + idx,   // 3
                    particles.v_ + idx        // 4
                    ));
        }

        else
        {
            return thrust::make_zip_iterator(       //
                thrust::make_tuple(                 //
                    particles.weight_.data() + idx, // 0
                    particles.charge_.data() + idx, // 1
                    particles.iCell_.data() + idx,  // 2
                    particles.delta_.data() + idx,  // 3
                    particles.v_.data() + idx       // 4
                    ));
        }
    }
    template<typename T>
    static auto& weight(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<0>(it);
    }
    template<typename T>
    static auto& charge(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<1>(it);
    }

    template<typename T>
    static auto& iCell(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<2>(it);
    }
    template<typename T>
    static auto& delta(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<3>(it);
    }
    template<typename T>
    static auto& v(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<4>(it);
    }
};


template<typename SoAParticles_t>
struct SoAZipParticle
{
    auto constexpr static dim = SoAParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAIteratorAdaptor::template make<SoAParticles_t>(
        std::declval<SoAParticles_t&>()))>;

    SoAZipParticle(SoAParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAIteratorAdaptor::make(ps, i)}
    {
    }

    auto& weight() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::weight(*it); }
    auto& weight() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::weight(*it); }
    auto& charge() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::charge(*it); }
    auto& charge() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::charge(*it); }
    auto& iCell() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::iCell(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::iCell(*it); }
    auto& delta() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::delta(*it); }
    auto& delta() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::delta(*it); }
    auto& v() _PHARE_ALL_FN_ { return SoAIteratorAdaptor::v(*it); }
    auto& v() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::v(*it); }

    Iterator it;
};



template<typename SoAParticles_t>
struct SoAZipConstParticle
{
    auto constexpr static dim = SoAParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAIteratorAdaptor::template make<SoAParticles_t>(
        std::declval<SoAParticles_t&>()))>;

    SoAZipConstParticle(SoAParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAIteratorAdaptor::make(ps, i)}
    {
    }

    auto& weight() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::weight(*it); }
    auto& charge() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::charge(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::iCell(*it); }
    auto& delta() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::delta(*it); }
    auto& v() const _PHARE_ALL_FN_ { return SoAIteratorAdaptor::v(*it); }

    auto copy() const _PHARE_ALL_FN_
    {
        return Particle<dim>{weight(), charge(), iCell(), delta(), v()};
    }

    Iterator it;
};


struct SoAVXIteratorAdaptor
{
    template<typename Particles>
    static auto soaxv_tuple_as_thrust_tuple(Particles& v, std::size_t const idx)
    {
        auto constexpr static dim = Particles::dimension;

        if constexpr (dim == 1)
            return thrust::make_tuple(&v.weight_[idx],   // 0
                                      &v.charge_[idx],   // 1
                                      &v.iCell_[idx],    // 2
                                      &v.delta_[0][idx], // 3
                                      &v.v_[0][idx],     // 4
                                      &v.v_[1][idx],     // 5
                                      &v.v_[2][idx]      // 6
            );

        if constexpr (dim == 2)
            return thrust::make_tuple(&v.weight_[idx],   // 0
                                      &v.charge_[idx],   // 1
                                      &v.iCell_[idx],    // 2
                                      &v.delta_[0][idx], // 3
                                      &v.delta_[1][idx], // 4
                                      &v.v_[0][idx],     // 5
                                      &v.v_[1][idx],     // 6
                                      &v.v_[2][idx]      // 7
            );

        if constexpr (dim == 3)
            return thrust::make_tuple(&v.weight_[idx],   // 0
                                      &v.charge_[idx],   // 1
                                      &v.iCell_[idx],    // 2
                                      &v.delta_[0][idx], // 3
                                      &v.delta_[1][idx], // 4
                                      &v.delta_[2][idx], // 5
                                      &v.v_[0][idx],     // 6
                                      &v.v_[1][idx],     // 7
                                      &v.v_[2][idx]      // 8
            );
    }

    template<typename Particles>
    auto static make(Particles& particles, std::size_t const idx = 0) _PHARE_ALL_FN_
    {
        auto static constexpr dim = Particles::dimension;

        return thrust::make_zip_iterator(soaxv_tuple_as_thrust_tuple(particles, idx));
    }


    template<typename T>
    static auto& weight(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<0>(it);
    }
    template<typename T>
    static auto& charge(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<1>(it);
    }

    template<typename T>
    static auto& iCell(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<2>(it);
    }
    template<std::uint8_t dim, std::uint8_t idx = 0, typename T>
    static auto& delta(T&& it) _PHARE_ALL_FN_
    {
        return thrust::get<3 + idx>(it);
    }

    template<std::uint8_t dim, typename T>
    static auto fdelta(T&& it) _PHARE_ALL_FN_
    {
        return for_N<dim, for_N_R_mode::make_array>(
            [&](auto d) { return SoAVXIteratorAdaptor::template delta<dim, d>(it); });
    }

    template<std::uint8_t dim, std::uint8_t idx = 0, typename T>
    static auto& v(T&& it) _PHARE_ALL_FN_
    {
        // 1d start = 4  // 3 + 1 + 0 == 4
        // 3d start = 6  // 3 + 3 + 0 == 6

        return thrust::get<3 + dim + idx>(it);
    }

    template<std::uint8_t dim, typename T>
    static auto fv(T&& it) _PHARE_ALL_FN_
    {
        return for_N<3, for_N_R_mode::make_array>(
            [&](auto d) { return SoAVXIteratorAdaptor::template v<dim, d>(it); });
    }
};


template<typename SoAVXParticles_t>
struct SoAVXZipParticle
{
    auto constexpr static dim = SoAVXParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAVXIteratorAdaptor::template make<SoAVXParticles_t>(
        std::declval<SoAVXParticles_t&>()))>;

    SoAVXZipParticle(SoAVXParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAVXIteratorAdaptor::make(ps, i)}
    {
    }

    void set(auto that) { it = that.it; }

    auto& weight() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::weight(*it); }
    auto& weight() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::weight(*it); }
    auto& charge() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::charge(*it); }
    auto& charge() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::charge(*it); }
    auto& iCell() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::iCell(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::iCell(*it); }
    // auto& delta() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::delta(*it); }
    // auto& delta() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::delta(*it); }
    // auto& v() _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::v(*it); }
    // auto& v() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::v(*it); }
    auto fdelta() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::template fdelta<dim>(*it); }
    auto fv() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::template fv<dim>(*it); }

    Iterator it;
};



template<typename SoAVXParticles_t>
struct SoAVXZipConstParticle
{
    auto constexpr static dim = SoAVXParticles_t::dimension;
    using Iterator = std::decay_t<decltype(SoAVXIteratorAdaptor::template make<SoAVXParticles_t>(
        std::declval<SoAVXParticles_t&>()))>;

    SoAVXZipConstParticle(SoAVXParticles_t& ps, std::size_t const& i) _PHARE_ALL_FN_
        : it{SoAVXIteratorAdaptor::make(ps, i)}
    {
    }
    void set(auto that) { it = that.it; }

    auto& weight() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::weight(*it); }
    auto& charge() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::charge(*it); }
    auto& iCell() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::iCell(*it); }
    // auto& delta() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::delta(*it); }
    // auto& v() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::v(*it); }

    auto fdelta() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::template fdelta<dim>(*it); }
    auto fv() const _PHARE_ALL_FN_ { return SoAVXIteratorAdaptor::template fv<dim>(*it); }

    auto copy() const _PHARE_ALL_FN_
    {
        return Particle<dim>{weight(), charge(), iCell(),
                             for_N<dim, for_N_R_mode::make_array>([&](auto d) {
                                 return SoAVXIteratorAdaptor::template delta<dim, d>(*it);
                             }),
                             for_N<3, for_N_R_mode::make_array>([&](auto d) -> auto& {
                                 return SoAVXIteratorAdaptor::template v<dim, d>(*it);
                             })};
    }

    Iterator it;
};




template<typename SoAParticles_t, bool _is_const = false>
struct SoAZipParticle_t
{
    bool static constexpr is_const
        = _is_const || std::is_const_v<std::remove_reference_t<SoAParticles_t>>;

    using value_type = std::conditional_t<is_const, SoAZipConstParticle<SoAParticles_t const>,
                                          SoAZipParticle<SoAParticles_t>>;
};

template<typename SoAParticles_t, bool _is_const = false>
struct SoAVXZipParticle_t
{
    bool static constexpr is_const
        = _is_const || std::is_const_v<std::remove_reference_t<SoAParticles_t>>;

    using value_type = std::conditional_t<is_const, SoAVXZipConstParticle<SoAParticles_t const>,
                                          SoAVXZipParticle<SoAParticles_t>>;
};


template<typename Particles>
auto particle_zip_iterator(Particles& ps, std::size_t const i)
{
    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, SoA))
        return typename SoAZipParticle_t<Particles>::value_type{ps, i};
    else if constexpr (any_in(Particles::layout_mode, SoAVX))
        return typename SoAVXZipParticle_t<Particles>::value_type{ps, i};
    else
        throw std::runtime_error("no impl");
}

template<typename T, std::size_t dim>
auto partitionner(SoAIteratorAdaptor& begin, SoAIteratorAdaptor& end, Box<T, dim> const& box)
{
    return partitionner(begin, end, box, [&box](auto const& part) {
        return isIn(SoAIteratorAdaptor::iCell(part), box);
    });
}


#endif //  PHARE_HAVE_THRUST

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_THRUST_HPP */
