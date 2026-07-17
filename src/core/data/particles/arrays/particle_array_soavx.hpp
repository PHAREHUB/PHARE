#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOAVX_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOAVX_HPP

#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/arrays/particle_array_soa_thrust.hpp"

#include <tuple>
#include <type_traits>

namespace PHARE::core
{

template<std::size_t dim>
static auto soaxv_as_tuple(auto& v)
{
    auto tup = std::tuple_cat(
        std::forward_as_tuple(v.weight_, v.charge_, v.iCell_),
        for_N<dim, for_N_R_mode::forward_tuple>([&](auto i) -> auto& { return v.delta_[i]; }),
        for_N<3, for_N_R_mode::forward_tuple>([&](auto i) -> auto& { return v.v_[i]; }) //
    );

    if constexpr (dim == 1)
    {
        static_assert(std::tuple_size_v<decltype(tup)> == 7);
    }
    if constexpr (dim == 2)
    {
        static_assert(std::tuple_size_v<decltype(tup)> == 8);
    }
    if constexpr (dim == 3)
    {
        static_assert(std::tuple_size_v<decltype(tup)> == 9);
    }

    return tup;
}

template<std::size_t dim, std::size_t size_>
struct SoAVXArray
{
    // using This                         = SoAVXArray<dim, size_>;
    auto static constexpr storage_mode = StorageMode::ARRAY;
    auto static constexpr dimension    = dim;
    auto static constexpr is_vector    = false;

    auto constexpr static size() { return size_; }

    template<typename Particles_t>
    void assign(Particles_t const& src, std::size_t const idx, std::size_t const dst) _PHARE_ALL_FN_
    {
        assert(dst < size_);
        assert(idx < src.size());
        this->weight_[dst] = src.weight(idx);
        this->charge_[dst] = src.charge(idx);
        this->iCell_[dst]  = src.iCell(idx);
        for_N<dim>([&](auto vi) { this->delta_[vi][dst] = src.delta(idx)[vi]; });
        for_N<3>([&](auto vi) { this->v_[vi][dst] = src.v(idx)[vi]; });
    }

    template<typename Particle_t>
    void assign(Particle_t const& src, std::size_t const dst) _PHARE_ALL_FN_
    {
        assert(dst < size_);
        this->weight_[dst] = src.weight();
        this->charge_[dst] = src.charge();
        this->iCell_[dst]  = src.iCell();
        for_N<dim>([&](auto vi) { this->delta_[vi][dst] = src.delta()[vi]; });
        for_N<3>([&](auto vi) { this->v_[vi][dst] = src.v()[vi]; });
    }

#if PHARE_HAVE_THRUST
    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return SoAVXZipParticle(*this, s); }

    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return SoAVXZipConstParticle(*this, s);
    }
#endif

    auto as_tuple() _PHARE_ALL_FN_ { return soaxv_as_tuple<dim>(*this); }
    auto as_tuple() const _PHARE_ALL_FN_ { return soaxv_as_tuple<dim>(*this); }

    std::array<double, size_> weight_, charge_;
    std::array<std::array<int, size_>, dim> iCell_;
    std::array<std::array<double, size_>, dim> delta_;
    std::array<std::array<double, size_>, 3> v_;
};




// used when the memory is owned elsewhere, e.g. numpy arrays
template<std::size_t dim, auto alloc_mode_>
struct SoAVXSpan
{
    static_assert(std::is_same_v<decltype(alloc_mode_), PHARE::AllocatorMode>);

    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = StorageMode::SPAN;
    auto static constexpr dimension    = dim;

    using SIZE_T = default_span_size_t;

    template<typename T>
    using container_t = T*;

    SoAVXSpan() = default;

    template<typename ParticleArray>
    SoAVXSpan(ParticleArray&& array, std::size_t const& beg, std::size_t const& siz)
        : size_{siz}
    {
        reset(array);
    }

    template<typename ParticleArray>
    SoAVXSpan(ParticleArray&& array, std::size_t const& siz)
        : SoAVXSpan{array, 0, siz}
    {
    }

    template<typename ParticleArray, // SFINAE protection to only allow particle arrays
             typename = std::enable_if_t<ParticleArray::layout_mode == LayoutMode::SoAVX>>
    SoAVXSpan(ParticleArray& array)
        : SoAVXSpan{array, 0, array.size()}
    {
    }

    template<typename ParticleArray, // SFINAE protection to only allow particle arrays
             typename = std::enable_if_t<ParticleArray::layout_mode == LayoutMode::SoAVX>>
    SoAVXSpan(ParticleArray&& array)
        : SoAVXSpan{array, 0, array.size()}
    {
    }

    auto size() const _PHARE_ALL_FN_ { return size_; }
    void clear() _PHARE_ALL_FN_ { size_ = 0; }
    void resize(std::size_t const& s)
    {
        // PHARE_ASSERT(s <= size_); // can't be bigger
        size_ = s;
    }

    void pop_back() _PHARE_ALL_FN_ { --size_; }
    auto size_address() _PHARE_ALL_FN_ { return &size_; }

    auto as_tuple() _PHARE_ALL_FN_ { return soaxv_as_tuple<dim>(*this); }
    auto as_tuple() const _PHARE_ALL_FN_ { return soaxv_as_tuple<dim>(*this); }

    template<typename Particles_t>
    void reset(Particles_t& src)
    {
        this->weight_ = src.weight_.data();
        this->charge_ = src.charge_.data();
        this->iCell_  = src.iCell_.data();
        for_N<dim>([&](auto vi) { this->delta_[vi] = src.delta()[vi].data(); });
        for_N<3>([&](auto vi) { this->v_[vi] = src.v()[vi].data(); });
        size_ = src.size();
    }

    SIZE_T size_;
    container_t<double> weight_, charge_;
    container_t<std::array<int, dim>> iCell_;
    std::array<container_t<double>, dim> delta_;
    std::array<container_t<double>, 3> v_;
};



template<std::size_t dim, auto alloc_mode_>
struct SoAVXVector
{
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr dimension    = dim;
    using This                         = SoAVXVector<dim, alloc_mode>;

    template<typename Type>
    using container_t = typename Vector<Type, alloc_mode, 1>::vector_t;

    SoAVXVector() {}

    SoAVXVector(std::size_t size)
        : weight_(size)
        , charge_(size)
        , iCell_(size)
        , delta_(size)
        , v_(size)
    {
    }

    auto size() const _PHARE_ALL_FN_ { return weight_.size(); }

    void pop_back()
    {
        std::apply([](auto&... v) { ((v.pop_back()), ...); }, as_tuple());
    }


    auto as_tuple() { return soaxv_as_tuple<dim>(*this); }
    auto as_tuple() const { return soaxv_as_tuple<dim>(*this); }

    void clear()
    {
        std::apply([](auto&... container) { ((container.clear()), ...); }, as_tuple());
    }

    void resize(std::size_t const& size)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(std::apply(
                [&](auto&... container) {
                    ((mkn::gpu::resize(container, size, /*copy=*/true)), ...);
                },
                as_tuple()));
        }
        else
            std::apply([&](auto&... container) { ((container.resize(size)), ...); }, as_tuple());
    }


    container_t<double> weight_, charge_;

    // NO AVX FOR ICELL DUE TO thrust::tuple zip_iterator only supports 10 values
    container_t<std::array<int, dim>> iCell_;

    std::array<container_t<double>, dim> delta_;
    std::array<container_t<double>, 3> v_;

    template<typename V>
    static auto& get_vec(V& v)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(return mkn::gpu::as_super(v));
        }
        else
        {
            return v;
        }
    }
};


template<typename Super_>
class SoAVXParticles : public Super_
{
    template<typename T>
    struct iterator_impl;

public:
    using Super                        = Super_;
    using This                         = SoAVXParticles<Super>;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr layout_mode  = LayoutMode::SoAVX;
    auto static constexpr storage_mode = Super::storage_mode;

    using SIZE_T = default_span_size_t;
    // using Particle_t = SoAVXParticle_crt<dimension>;
    using Super::size;

    using Span_t = SoAVXParticles<SoAVXSpan<dimension, alloc_mode>>;
    friend class SoAVXParticles<SoAVXSpan<dimension, alloc_mode>>;

    // template<std::size_t size>
    //  using array_type = SoAVXParticles<SoAVXArray<dimension, size>>;

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    // public for pybind but avoid otherwise
    using Super::charge_;
    using Super::delta_;
    using Super::iCell_;
    using Super::weight_;

    using Super::v_;


    template<typename... Args>
    SoAVXParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    SoAVXParticles(iterator start, iterator end); // impl @ bottom of this file

    SoAVXParticles(This const& that)  = default;
    SoAVXParticles(This&& that)       = default;
    This& operator=(This&& that)      = default;
    This& operator=(This const& that) = default;


    auto& weight(std::size_t i) const _PHARE_ALL_FN_ { return weight_[i]; }
    auto& weight(std::size_t i) _PHARE_ALL_FN_ { return weight_[i]; }

    auto& charge(std::size_t i) const _PHARE_ALL_FN_ { return charge_[i]; }
    auto& charge(std::size_t i) _PHARE_ALL_FN_ { return charge_[i]; }

    auto& iCell(std::size_t i) const _PHARE_ALL_FN_ { return iCell_[i]; }
    auto& iCell(std::size_t i) _PHARE_ALL_FN_ { return iCell_[i]; }

    auto& delta(std::size_t i) const _PHARE_ALL_FN_ { return delta_[i]; }
    auto& delta(std::size_t i) _PHARE_ALL_FN_ { return delta_[i]; }

    auto& v(std::size_t i) const _PHARE_ALL_FN_ { return v_[i]; }
    auto& v(std::size_t i) _PHARE_ALL_FN_ { return v_[i]; }

    auto& weight() _PHARE_ALL_FN_ { return weight_; }
    auto& charge() _PHARE_ALL_FN_ { return charge_; }
    auto& iCell() _PHARE_ALL_FN_ { return iCell_; }
    auto& delta() _PHARE_ALL_FN_ { return delta_; }
    auto& v() _PHARE_ALL_FN_ { return v_; }

    auto& weight() const _PHARE_ALL_FN_ { return weight_; }
    auto& charge() const _PHARE_ALL_FN_ { return charge_; }
    auto& iCell() const _PHARE_ALL_FN_ { return iCell_; }
    auto& delta() const _PHARE_ALL_FN_ { return delta_; }
    auto& v() const _PHARE_ALL_FN_ { return v_; }


    // for performing the same operation across all vectors e.g. with std apply
    // auto as_tuple(std::size_t i) _PHARE_ALL_FN_
    // {
    //     return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
    //                                  this->delta_[i], this->v_[i]);
    // }

    // auto as_tuple(std::size_t i) const _PHARE_ALL_FN_
    // {
    //     return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
    //                                  this->delta_[i], this->v_[i]);
    // }

    auto as_tuple() _PHARE_ALL_FN_ { return soaxv_as_tuple<dimension>(*this); }
    auto as_tuple() const _PHARE_ALL_FN_ { return soaxv_as_tuple<dimension>(*this); }


    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args) _PHARE_ALL_FN_
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T*>{t, args...};
    }

    auto begin() const _PHARE_ALL_FN_ { return it(this); }
    auto begin() _PHARE_ALL_FN_ { return it(this); }
    auto end() const _PHARE_ALL_FN_ { return it(this, size()); }
    auto end() _PHARE_ALL_FN_ { return it(this, size()); }

    template<typename Impl, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void push_back(iterator_impl<Impl> const& iter)
    {
        Super::get_vec(this->weight_).push_back(iter.weight());
        Super::get_vec(this->charge_).push_back(iter.charge());
        Super::get_vec(this->iCell_).push_back(iter.iCell());
        for_N<dimension, for_N_R_mode::forward_tuple>([&](auto i) {
            Super::get_vec(this->delta_[i]).push_back(iter.template delta<dimension, i>());
        });
        for_N<3, for_N_R_mode::forward_tuple>([&](auto i) {
            Super::get_vec(this->v_[i]).push_back(iter.template v<dimension, i>());
        });
    }


    template<typename Particle_t, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void push_back(Particle_t const& particle)
    {
        Super::get_vec(this->weight_).push_back(particle.weight());
        Super::get_vec(this->charge_).push_back(particle.charge());
        Super::get_vec(this->iCell_).push_back(particle.iCell());
        for_N<dimension, for_N_R_mode::forward_tuple>(
            [&](auto i) { Super::get_vec(this->delta_[i]).push_back(particle.delta()[i]); });
        for_N<3, for_N_R_mode::forward_tuple>(
            [&](auto i) { Super::get_vec(this->v_[i]).push_back(particle.v()[i]); });
    }


    //     template<typename Impl, auto S = storage_mode,
    //              typename = std::enable_if_t<S == StorageMode::VECTOR>>
    //     void push_back(iterator_impl<Impl> const& particle)
    //     {
    //         Super::get_vec(this->weight_).push_back(particle.weight());
    //         Super::get_vec(this->charge_).push_back(particle.charge());
    //         Super::get_vec(this->iCell_).push_back(particle.iCell());
    //         Super::get_vec(this->delta_).push_back(particle.delta());
    //         Super::get_vec(this->v_).push_back(particle.v());
    //     }


    //     template<typename Particle_t, auto S = storage_mode,
    //              typename = std::enable_if_t<S == StorageMode::VECTOR>>
    //     void push_back(Particle_t const& particle)
    //     {
    //         auto const& [w, c, i, d, v] = particle;
    //         Super::get_vec(this->weight_).push_back(w);
    //         Super::get_vec(this->charge_).push_back(c);
    //         Super::get_vec(this->iCell_).push_back(i);
    //         Super::get_vec(this->delta_).push_back(d);
    //         Super::get_vec(this->v_).push_back(v);
    //     }


    // #if PHARE_HAVE_THRUST // needs thrust or no compile on non-const access

    //     template<typename Particle_t>
    //     auto emplace_back_zip(Particle_t const& particle)
    //     {
    //         Super::get_vec(this->weight_).emplace_back(particle.weight());
    //         Super::get_vec(this->charge_).emplace_back(particle.charge());
    //         Super::get_vec(this->iCell_).emplace_back(particle.iCell());
    //         Super::get_vec(this->delta_).emplace_back(particle.delta());
    //         Super::get_vec(this->v_).emplace_back(particle.v());
    //         return end() - 1;
    //     }

    //     template<typename That>
    //     auto emplace_back(SoAVXZipParticle<That> const& particle)
    //     {
    //         return emplace_back_zip(particle);
    //     }
    //     template<typename That>
    //     auto emplace_back(SoAVXZipConstParticle<That> const& particle)
    //     {
    //         return emplace_back_zip(particle);
    //     }

    // #endif


    //     template<typename Particle_t, auto S = storage_mode,
    //              typename = std::enable_if_t<S == StorageMode::VECTOR>>
    //     auto emplace_back(Particle_t const& particle)
    //     {
    //         Super::get_vec(this->weight_).emplace_back(particle.weight());
    //         Super::get_vec(this->charge_).emplace_back(particle.charge());
    //         Super::get_vec(this->iCell_).emplace_back(particle.iCell());
    //         Super::get_vec(this->delta_).emplace_back(particle.delta());
    //         Super::get_vec(this->v_).emplace_back(particle.v());
    //         return end() - 1;
    //     }

    //     template<typename That, auto S = storage_mode,
    //              typename = std::enable_if_t<S == StorageMode::VECTOR
    //                                          and That::layout_mode == LayoutMode::SoAVX>>
    //     void emplace_back(That const& src)
    //     {
    //         auto this_tuple = as_tuple();
    //         auto that_tuple = src.as_tuple();
    //         for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto vi) {
    //             auto& vec            = Super::get_vec(std::get<vi>(this_tuple));
    //             auto const& that_vec = std::get<vi>(that_tuple);
    //             for (std::size_t i = 0; i < src.size(); ++i)
    //                 vec.emplace_back(that_vec[i]);
    //         });
    //     }


    //     template<typename That, auto S = storage_mode,
    //              typename = std::enable_if_t<S == StorageMode::VECTOR
    //                                          and That::layout_mode == LayoutMode::SoAVX>>
    //     void emplace_back(That const& src, std::size_t const& idx)
    //     {
    //         auto this_tuple = as_tuple();
    //         auto that_tuple = src.as_tuple();
    //         for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto vi) {
    //             Super::get_vec(std::get<vi>(this_tuple)).emplace_back(std::get<vi>(that_tuple)[idx]);
    //         });
    //     }

    //     template<typename... Args>
    //     auto emplace_back(Args const&... args)
    //     {
    //         auto arg_tuple = std::forward_as_tuple(args...);

    //         if constexpr (std::tuple_size_v<decltype(arg_tuple)> == 5)
    //         {
    //             auto this_tuple = as_tuple();
    //             for_N<std::tuple_size_v<decltype(arg_tuple)>>([&](auto ic) {
    //                 auto constexpr i = ic();
    //                 Super::get_vec(std::get<i>(this_tuple)).emplace_back(std::get<i>(arg_tuple));
    //             });
    //         }

    //         else
    //             throw std::runtime_error("NO IMPL!");

    //         return end() - 1;
    //     }



    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void reserve(std::size_t const& size)
    {
        if constexpr (CompileOptions::WithMknGpu
                      and Super::alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(std::apply(
                [&](auto&... container) {
                    ((mkn::gpu::reserve(container, size, /*copy=*/true)), ...);
                },
                as_tuple()));
        }
        else
            std::apply([&](auto&... container) { ((container.reserve(size)), ...); }, as_tuple());
    }

    template<typename Src>
    void append(Src const& src, std::size_t const start, std::size_t const size)
    {
        reserve(this->size() + size);
        for (std::size_t i = 0; i < size; ++i)
            emplace_back(src[i + start]);
    }


    void swap(This& that)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = that.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto i) { std::get<i>(this_tuple).swap(std::get<i>(that_tuple)); });
    }

    void swap(std::size_t const& a, std::size_t const& b)
    {
        if (a == b)
            return;

        std::swap(weight_[a], weight_[b]);
        std::swap(charge_[a], charge_[b]);
        std::swap(iCell_[a], iCell_[b]);
        std::swap(delta_[a], delta_[b]);
        std::swap(v_[a], v_[b]);
    }




    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void erase(iterator_impl<This*> first, iterator_impl<This*> last)
    {
        std::apply(
            [&](auto&... container) {
                ((container.erase(container.begin() + first.curr_pos,
                                  container.begin() + last.curr_pos)),
                 ...);
            },
            as_tuple());
    }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    std::size_t capacity() const
    {
        return weight_.capacity(); // they're all the same
    }


    auto copy(std::size_t i) const _PHARE_ALL_FN_
    {
        return Particle<dimension>{
            weight_[i], charge_[i], iCell_[i],
            for_N<dimension, for_N_R_mode::make_array>(
                [&](auto d) -> auto& { return delta_[d][i]; }),
            for_N<3, for_N_R_mode::make_array>([&](auto d) -> auto& { return v_[d][i]; })};
    }


    auto back() { return (*this)[size() - 1]; }
    auto front() { return (*this)[0]; }


    auto constexpr static size_of_particle()
    {
        return sizeof(typename decltype(iCell_)::value_type)
               + sizeof(typename decltype(delta_)::value_type)
               + sizeof(typename decltype(weight_)::value_type)
               + sizeof(typename decltype(charge_)::value_type)
               + sizeof(typename decltype(v_)::value_type);
    }


    void check() const {}

    void assign(std::size_t const& src, std::size_t const& dst) _PHARE_ALL_FN_
    {
        std::apply([&](auto&... v) { ((v[dst] = v[src]), ...); }, as_tuple());
    }
    void assign(SIZE_T const& src, std::size_t const& dst) _PHARE_ALL_FN_
    {
        std::apply([&](auto&... v) { ((v[dst] = v[src]), ...); }, as_tuple());
    }

    void assign(Particle<dimension> const& src, std::size_t const& dst) _PHARE_ALL_FN_
    {
        auto this_tuple = as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto i) { std::get<i>(this_tuple)[dst] = std::get<i>(*src); });
    }

    template<typename Particle_t>
    void assign(Particle_t const& src, std::size_t const& dst) _PHARE_ALL_FN_
    {
        auto this_tuple = as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto i) { std::get<i>(this_tuple)[dst] = std::get<i>(*src); });
    }

    template<typename _Particles>
    void assign(_Particles const& src, std::size_t const& idx,
                std::size_t const& dst) _PHARE_ALL_FN_
    {
        auto this_tuple = as_tuple();
        auto that_tuple = src.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto vi) { std::get<vi>(this_tuple)[dst] = std::get<vi>(that_tuple)[idx]; });
    }

#if PHARE_HAVE_THRUST // needs thrust or no compile on non-const access
    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return SoAVXZipParticle(*this, s); }
#endif

    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return PHARE_WITH_THRUST(SoAVXZipConstParticle(*this, s));
        // else
        PHARE_WITH_THRUST_ELSE(copy(s));
    }
};


template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using SoAVXVectorParticles = SoAVXParticles<SoAVXVector<dim, alloc_mode>>;

template<std::size_t dim, std::size_t size>
using SoAVXArrayParticles = SoAVXParticles<SoAVXArray<dim, size>>;


#if PHARE_HAVE_THRUST
template<typename OuterSuper>
template<typename T>
struct SoAVXParticles<OuterSuper>::iterator_impl
{
    auto static constexpr dimension = OuterSuper::dimension;
    auto static constexpr is_const  = std::is_const_v<T>;
    auto static constexpr is_soavx  = true; // used for SFINAE

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using Particle_t =
        typename SoAVXZipParticle_t<SoAVXParticles<OuterSuper>, is_const>::value_type;
    using value_type = Particle_t;
    using pointer    = Particle_t*;
    using reference  = Particle_t&;


    iterator_impl(T& particles_, std::size_t const& s = 0) _PHARE_ALL_FN_
        : particles{particles_},
          curr_pos{s},
          particle{deref(particles), curr_pos}
    {
    }
    iterator_impl(iterator_impl&& that)      = default;
    iterator_impl(iterator_impl const& that) = default;

    auto& set()
    {
        // PHARE_LOG_LINE_SS(particle.fdelta()[0]);
        if (curr_pos != deref(particles).size())
            particle.set(particle_zip_iterator(deref(particles), curr_pos));
        // PHARE_LOG_LINE_SS(particle.fdelta()[0]);
        return *this;
    }

    auto& operator++()
    {
        ++curr_pos;
        return set();
    }
    auto operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    auto& operator+=(std::int64_t i) _PHARE_ALL_FN_
    {
        curr_pos += i;
        return set();
    }


    auto& operator--()
    {
        --curr_pos;
        return set();
    }

    // auto operator+(std::size_t i) const;

    auto operator+(std::int64_t i) const _PHARE_ALL_FN_
    {
        auto copy = *this;
        copy += i;
        return copy;
    }
    auto operator-(std::int64_t i) const
    {
        auto copy = *this;
        copy -= i;
        return copy;
    }

    auto& operator=(iterator_impl const& that)
    {
        curr_pos = that.curr_pos;
        return this->set();
    }
    auto& operator=(iterator_impl&& that)
    {
        curr_pos = that.curr_pos;
        return this->set();
    }

    auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
    auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }

    // auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
    // auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }
    // auto operator<(iterator_impl const& that) const { return curr_pos < that.curr_pos; }
    // auto operator-(iterator_impl const& that) { return curr_pos - that.curr_pos; }
    // auto& operator()() { return deref(particles); }
    // auto& operator()() const { return deref(particles); }
    auto& operator*() _PHARE_ALL_FN_ { return *this; }
    auto& operator*() const _PHARE_ALL_FN_ { return *this; }

    // auto idx() const { return curr_pos; }



    auto copy() const _PHARE_ALL_FN_ { return deref(particles).copy(curr_pos); }

    T particles;
    std::size_t curr_pos = 0;

    auto& weight() _PHARE_ALL_FN_ { return particle.weight(); }
    auto& weight() const _PHARE_ALL_FN_ { return particle.weight(); }
    auto& charge() _PHARE_ALL_FN_ { return particle.charge(); }
    auto& charge() const _PHARE_ALL_FN_ { return particle.charge(); }
    // auto& iCell() _PHARE_ALL_FN_ { return deref(particles).iCell_[curr_pos]; }

    auto& iCell() _PHARE_ALL_FN_ { return particle.iCell(); }
    auto& iCell() const _PHARE_ALL_FN_ { return particle.iCell(); }

    // auto& delta() _PHARE_ALL_FN_ { return deref(particles).delta_[curr_pos]; }
    auto delta() const _PHARE_ALL_FN_ { return particle.fdelta(); }
    // auto& v() _PHARE_ALL_FN_ { return deref(particles).v_[curr_pos]; }
    auto v() const _PHARE_ALL_FN_ { return particle.fv(); }

    Particle_t particle;
    // auto particle() const { return particle_zip_iterator(deref(particles), curr_pos); }
};
#endif // PHARE_HAVE_THRUST



} // namespace PHARE::core


namespace std
{

// template<template<typename> typename Iterator, typename O>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::is_soa, PHARE::core::Particle<Iterator<O>::dimension>>
// copy(Iterator<O> src)
// {
//     return {src.weight(), src.charge(), src.iCell(), src.delta(), src.v()};
// }


// template<template<typename> typename Iterator, typename O>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::is_soa, void> copy(Iterator<O> src_begin, //
//                                                           Iterator<O> src_end,
//                                                           Iterator<O> dst_begin)
// {
//     auto src_tuple = src_begin.particles.as_tuple();
//     auto dst_tuple = dst_begin.particles.as_tuple();
//     for_N<std::tuple_size_v<decltype(src_tuple)>>([&](auto i) {
//         auto& src = std::get<i>(src_tuple);
//         auto& dst = std::get<i>(dst_tuple);
//         std::copy(&src[src_begin.curr_pos], &src[src_end.curr_pos],
//         &dst[dst_begin.curr_pos]);
//     });
// }



// template<template<typename> typename Iterator, typename O, typename Inserter>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::is_soa, void> copy(Iterator<O> src_begin, //
//                                                           Iterator<O> src_end, Inserter i)
// {
//     // auto const& particles = src_begin.particles;

//     for (; src_begin != src_end; ++src_begin, ++i)
//         *i = src_begin.copy();
// }

// template<typename OuterSuper, typename T, typename Inserter, typename Fn,
//          typename Iter
//          = typename PHARE::core::SoAVXParticles<OuterSuper>::template iterator_impl<T>>
// void copy_if(Iter src_begin, Iter src_end, Inserter i, Fn check)
// {
//     auto const& particles = src_begin.particles;

//     for (; src_begin != src_end; ++src_begin)
//         if (check(src_begin))
//         {
//             *i = src_begin.copy();
//             ++i;
//         }
// }



// template<typename OuterSuper, typename T,
//          typename Iter
//          = typename PHARE::core::SoAVXParticles<OuterSuper>::template iterator_impl<T>>
// void iter_swap(Iter a, Iter b)
// {
//     PHARE_LOG_LINE_STR("");
//     PHARE_ASSERT(&a.particles == &b.particles);

//     a.particles.swap(a.curr_pos, b.curr_pos);
// }

// template<typename OuterSuper, typename T, typename Predicate,
//          typename Iter
//          = typename PHARE::core::SoAVXParticles<OuterSuper>::template iterator_impl<T>>
// Iter find_if_not(Iter first, Iter last, Predicate q)
// { // https://en.cppreference.com/w/cpp/algorithm/find
//     PHARE_LOG_LINE_STR("");
//     for (; first != last; ++first)
//         if (!q(first)) // pass iterator!
//             return first;
//     return last;
// }

// template<typename OuterSuper, typename T, typename Predicate/*,
//          typename Iter = typename PHARE::core::SoAVXParticles<OuterSuper>::template
//          iterator_impl<T>*/>
// typename PHARE::core::SoAVXParticles<OuterSuper>::template iterator_impl<T>
// partition(typename PHARE::core::SoAVXParticles<OuterSuper>::template iterator_impl<T> first,
// typename PHARE::core::SoAVXParticles<OuterSuper>::template iterator_impl<T> last, Predicate
// p) { // https://en.cppreference.com/w/cpp/algorithm/partition
//     PHARE_LOG_LINE_STR("");
//     first = find_if_not(first, last, p);
//     if (first == last)
//         return first;

//     for (auto i = first + 1; i != last; ++i)
//     {
//         if (p(i)) // pass iterator!
//         {
//             iter_swap(i, first);
//             ++first;
//         }
//     }
//     return first;
// }




// template<template<typename> typename Iterator, typename O, typename Fn>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::is_soa, void>
// transform(Iterator<O> src_begin, Iterator<O> src_end, Iterator<O> dst_begin, Fn transform)
// {
//     if (&src_begin.particles != &dst_begin.particles)
//         throw std::runtime_error("transform cannot work");

//     for (; src_begin != src_end; ++src_begin, ++dst_begin)
//         /*dst_begin = */ transform(src_begin); // transform edits in place
// }


// template<template<typename> typename Iterator, typename O>
// // SFINAE to support const iterators
// typename std::enable_if_t<Iterator<O>::is_soavx, std::size_t> distance(Iterator<O> const& a,
//                                                                        Iterator<O> const& b)
// {
//     PHARE_ASSERT(&deref(a.particles).weight(0) == &deref(b.particles).weight(0));
//     PHARE_ASSERT(a.curr_pos <= b.curr_pos);

//     return b.curr_pos - a.curr_pos;
// }

} // namespace std


namespace PHARE::core
{


template<typename Super>
SoAVXParticles<Super>::SoAVXParticles(SoAVXParticles<Super>::iterator start,
                                      SoAVXParticles<Super>::iterator end)
    : Super{std::distance(start, end)}
{
    std::copy(start, end, this->begin()); // impl above
}

} // namespace PHARE::core




#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOAVX_HPP */
