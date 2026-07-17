#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP

#include "core/def/detail/mkn_avx.hpp"
#include "core/logger.hpp"
#include "core/vector.hpp"
// #include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "core/data/particles/arrays/particle_array_soa_thrust.hpp"


namespace PHARE::core
{


template<std::size_t dim, std::size_t size_>
struct SoAArray
{
    auto static constexpr storage_mode = StorageMode::ARRAY;
    auto static constexpr dimension    = dim;
    auto static constexpr is_vector    = false;

    std::array<double, size_> weight_, charge_;
    std::array<std::array<int, dim>, size_> iCell_;
    std::array<std::array<double, dim>, size_> delta_;
    std::array<std::array<double, 3>, size_> v_;

    auto constexpr static size() { return size_; }

    template<typename Particles_t>
    void assign(Particles_t const& src, std::size_t const idx, std::size_t const dst) _PHARE_ALL_FN_
    {
        assert(dst < size_);
        assert(idx < src.size());
        this->weight_[dst] = src.weight(idx);
        this->charge_[dst] = src.charge(idx);
        this->iCell_[dst]  = src.iCell(idx);
        this->delta_[dst]  = src.delta(idx);
        this->v_[dst]      = src.v(idx);
    }

    template<typename Particle_t>
    void assign(Particle_t const& src, std::size_t const dst) _PHARE_ALL_FN_
    {
        assert(dst < size_);
        this->weight_[dst] = src.weight();
        this->charge_[dst] = src.charge();
        this->iCell_[dst]  = src.iCell();
        this->delta_[dst]  = src.delta();
        this->v_[dst]      = src.v();
    }

#if PHARE_HAVE_THRUST

    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return SoAZipParticle(*this, s); }

    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return SoAZipConstParticle(*this, s);
    }

#endif

    auto as_tuple() _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }
    auto as_tuple() const _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }
};



// used when the memory is owned elsewhere, e.g. numpy arrays
template<std::size_t dim, auto alloc_mode_, auto _const_ = 0>
struct SoASpan
{
    static_assert(std::is_same_v<decltype(alloc_mode_), PHARE::AllocatorMode>);

    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = StorageMode::SPAN;
    auto static constexpr dimension    = dim;
    template<typename T>
    using container_t = std::conditional_t<_const_, T const*, T*>;
    using SIZE_T      = default_span_size_t;

    SoASpan() = default;

    // template<typename C0, typename C1>
    SoASpan(auto&& _iCell, auto&& _delta, auto&& _weight, auto&& _charge, auto&& _v)
        : size_{_weight.size()}
        , weight_{reinterpret_cast<container_t<double>>(_weight.data())}
        , charge_{reinterpret_cast<container_t<double>>(_charge.data())}
        , iCell_{reinterpret_cast<container_t<std::array<int, dim>>>(_iCell.data())}
        , delta_{reinterpret_cast<container_t<std::array<double, dim>>>(_delta.data())}
        , v_{reinterpret_cast<container_t<std::array<double, 3>>>(_v.data())}
    {
    }


    template<typename ParticleArray>
    SoASpan(ParticleArray&& array, std::size_t const& beg, std::size_t const& siz)
        : size_{siz}
        , weight_{&array.weight_[0] + beg}
        , charge_{&array.charge_[0] + beg}
        , iCell_{&array.iCell_[0] + beg}
        , delta_{&array.delta_[0] + beg}
        , v_{&array.v_[0] + beg}
    {
    }

    template<typename ParticleArray>
    SoASpan(ParticleArray&& array, std::size_t const& siz)
        : SoASpan{array, 0, siz}
    {
    }

    template<typename ParticleArray, // SFINAE protection to only allow particle arrays
             typename = std::enable_if_t<ParticleArray::layout_mode == LayoutMode::SoA>>
    SoASpan(ParticleArray&& array)
        : SoASpan{array, 0, array.size()}
    {
    }
    template<typename ParticleArray, // SFINAE protection to only allow particle arrays
             typename = std::enable_if_t<ParticleArray::layout_mode == LayoutMode::SoA>>
    SoASpan(ParticleArray& array)
        : SoASpan{array, 0, array.size()}
    {
    }

    auto size() const _PHARE_ALL_FN_ { return size_; }
    void clear() _PHARE_ALL_FN_ { size_ = 0; }
    void resize(std::size_t const& s)
    {
        PHARE_ASSERT(s <= size_); // can't be bigger
        size_ = s;
    }


    void pop_back() _PHARE_ALL_FN_ { --size_; }
    auto size_address() _PHARE_ALL_FN_ { return &size_; }


    auto as_tuple() _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }
    auto as_tuple() const _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }

    template<typename Particles_t>
    void reset(Particles_t& particles)
    {
        auto self = as_tuple();
        auto that = particles.as_tuple();
        for_N<std::tuple_size_v<decltype(that)>>(
            [&](auto vi) { std::get<vi>(self) = std::get<vi>(that).data(); });
        size_ = particles.size();
    }

    SIZE_T size_;
    container_t<double> weight_, charge_;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta_;
    container_t<std::array<double, 3>> v_;
};



template<std::size_t dim, auto alloc_mode_>
struct SoAVector
{
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr dimension    = dim;

    std::uint8_t constexpr static alloc_impl()
    {
        if (alloc_mode_ == AllocatorMode::GPU_UNIFIED)
            return 1;
        return PHARE_HAVE_MKN_AVX;
    }


    template<typename Type>
    using container_t = typename Vector<Type, alloc_mode, alloc_impl()>::vector_t;

    SoAVector() {}

    SoAVector(std::size_t size)
        : weight_(size)
        , charge_(size)
        , iCell_(size)
        , delta_(size)
        , v_(size)
    {
    }

    template<typename Particle_t>
    SoAVector(std::size_t size, Particle_t&& from)
        : weight_(size, from.weight())
        , charge_(size, from.charge())
        , iCell_(size, from.iCell())
        , delta_(size, from.delta())
        , v_(size, from.v())
    {
    }

    auto size() const _PHARE_ALL_FN_ { return weight_.size(); }


    void pop_back()
    {
        std::apply([](auto&... v) { ((v.pop_back()), ...); }, as_tuple());
    }

    auto as_tuple() { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
    auto as_tuple() const { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }

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
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta_;
    container_t<std::array<double, 3>> v_;

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
class SoAParticles : public Super_
{
    template<typename T>
    struct iterator_impl;

public:
    using Super                        = Super_;
    using This                         = SoAParticles<Super>;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr layout_mode  = LayoutMode::SoA;
    auto static constexpr storage_mode = Super::storage_mode;

    using SIZE_T     = default_span_size_t;
    using Particle_t = SoAParticle_crt<dimension>;
    using Super::size;

    using Span_t = SoAParticles<SoASpan<dimension, alloc_mode>>;
    friend class SoAParticles<SoASpan<dimension, alloc_mode>>;

    // template<std::size_t size>
    //  using array_type = SoAParticles<SoAArray<dimension, size>>;

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    // public for pybind but avoid otherwise
    using Super::charge_;
    using Super::delta_;
    using Super::iCell_;
    // using Super::resize;
    using Super::v_;
    using Super::weight_;


    template<typename... Args>
    SoAParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    SoAParticles(iterator start, iterator end); // impl @ bottom of this file

    SoAParticles(This const& that)    = default;
    SoAParticles(This&& that)         = default;
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
    auto as_tuple(std::size_t i) _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i]);
    }

    auto as_tuple(std::size_t i) const _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i]);
    }

    auto as_tuple() _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }
    auto as_tuple() const _PHARE_ALL_FN_
    {
        return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }


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
    void push_back(iterator_impl<Impl> const& particle)
    {
        Super::get_vec(this->weight_).push_back(particle.weight());
        Super::get_vec(this->charge_).push_back(particle.charge());
        Super::get_vec(this->iCell_).push_back(particle.iCell());
        Super::get_vec(this->delta_).push_back(particle.delta());
        Super::get_vec(this->v_).push_back(particle.v());
    }


    template<typename Particle_t, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void push_back(Particle_t const& particle)
    {
        auto const& [w, c, i, d, v] = particle;
        Super::get_vec(this->weight_).push_back(w);
        Super::get_vec(this->charge_).push_back(c);
        Super::get_vec(this->iCell_).push_back(i);
        Super::get_vec(this->delta_).push_back(d);
        Super::get_vec(this->v_).push_back(v);
    }


#if PHARE_HAVE_THRUST // needs thrust or no compile on non-const access

    template<typename Particle_t>
    auto emplace_back_zip(Particle_t const& particle)
    {
        Super::get_vec(this->weight_).emplace_back(particle.weight());
        Super::get_vec(this->charge_).emplace_back(particle.charge());
        Super::get_vec(this->iCell_).emplace_back(particle.iCell());
        Super::get_vec(this->delta_).emplace_back(particle.delta());
        Super::get_vec(this->v_).emplace_back(particle.v());
        return end() - 1;
    }

    template<typename That>
    auto emplace_back(SoAZipParticle<That> const& particle)
    {
        return emplace_back_zip(particle);
    }
    template<typename That>
    auto emplace_back(SoAZipConstParticle<That> const& particle)
    {
        return emplace_back_zip(particle);
    }

#endif


    template<typename Particle_t, auto S = storage_mode,
             typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto emplace_back(Particle_t const& particle)
    {
        Super::get_vec(this->weight_).emplace_back(particle.weight());
        Super::get_vec(this->charge_).emplace_back(particle.charge());
        Super::get_vec(this->iCell_).emplace_back(particle.iCell());
        Super::get_vec(this->delta_).emplace_back(particle.delta());
        Super::get_vec(this->v_).emplace_back(particle.v());
        return end() - 1;
    }

    template<typename That, auto S = storage_mode,
             typename
             = std::enable_if_t<S == StorageMode::VECTOR and That::layout_mode == LayoutMode::SoA>>
    void emplace_back(That const& src)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = src.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto vi) {
            auto& vec            = Super::get_vec(std::get<vi>(this_tuple));
            auto const& that_vec = std::get<vi>(that_tuple);
            for (std::size_t i = 0; i < src.size(); ++i)
                vec.emplace_back(that_vec[i]);
        });
    }


    template<typename That, auto S = storage_mode,
             typename
             = std::enable_if_t<S == StorageMode::VECTOR and That::layout_mode == LayoutMode::SoA>>
    void emplace_back(That const& src, std::size_t const& idx)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = src.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto vi) {
            Super::get_vec(std::get<vi>(this_tuple)).emplace_back(std::get<vi>(that_tuple)[idx]);
        });
    }

    template<typename... Args>
    auto emplace_back(Args const&... args)
    {
        auto arg_tuple = std::forward_as_tuple(args...);

        if constexpr (std::tuple_size_v<decltype(arg_tuple)> == 5)
        {
            auto this_tuple = as_tuple();
            for_N<std::tuple_size_v<decltype(arg_tuple)>>([&](auto ic) {
                auto constexpr i = ic();
                Super::get_vec(std::get<i>(this_tuple)).emplace_back(std::get<i>(arg_tuple));
            });
        }

        else
            throw std::runtime_error("NO IMPL!");

        return end() - 1;
    }



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

    // template<typename Src>
    // void append(Src const& src, std::size_t const start, std::size_t const size)
    // {
    //     if (capacity() < this->size() + size)
    //         reserve(this->size() + size);
    //     for (std::size_t i = 0; i < size; ++i)
    //         emplace_back(src[i + start]);
    // }

    // template<typename Src>
    // void append(Src const& src)
    // {
    //     append(src, 0, src.size());
    // }


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
            weight_[i], charge_[i], iCell_[i], delta_[i], v_[i],
        };
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
    auto operator[](std::size_t const& s) _PHARE_ALL_FN_ { return SoAZipParticle(*this, s); }
#endif

    auto operator[](std::size_t const& s) const _PHARE_ALL_FN_
    {
        return PHARE_WITH_THRUST(SoAZipConstParticle(*this, s));
        // else
        PHARE_WITH_THRUST_ELSE(copy(s));
    }
};


template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using SoAVectorParticles = SoAParticles<SoAVector<dim, alloc_mode>>;

template<std::size_t dim, std::size_t size>
using SoAArrayParticles = SoAParticles<SoAArray<dim, size>>;



template<typename OuterSuper>
template<typename T>
struct SoAParticles<OuterSuper>::iterator_impl
{
    auto static constexpr dimension = OuterSuper::dimension;
    // auto static constexpr is_const  = std::is_const_v<T>;
    auto static constexpr is_soa = true; // used for SFINAE

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using Particle_t        = typename ParticleDefaults<dimension>::Particle_t;
    using value_type        = Particle_t;
    using pointer           = Particle_t*;
    using reference         = Particle_t&;

    iterator_impl(T& particles_, std::size_t const& s = 0) _PHARE_ALL_FN_ : particles{particles_},
                                                                            curr_pos{s}
    {
    }
    iterator_impl(iterator_impl&& that)      = default;
    iterator_impl(iterator_impl const& that) = default;

    auto& operator++()
    {
        ++curr_pos;
        return *this;
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
        return *this;
    }


    auto& operator--()
    {
        --curr_pos;
        return *this;
    }
    auto operator+(std::int64_t i) const _PHARE_ALL_FN_
    {
        auto copy = *this;
        copy.curr_pos += i;
        return copy;
    }
    auto operator-(std::int64_t i) const
    {
        auto copy = *this;
        copy.curr_pos -= i;
        return copy;
    }

    auto& operator=(iterator_impl const& that)
    {
        curr_pos = that.curr_pos;
        return *this;
    }
    auto& operator=(iterator_impl&& that)
    {
        curr_pos = that.curr_pos;
        return *this;
    }


    auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
    auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }
    auto operator<(iterator_impl const& that) const { return curr_pos < that.curr_pos; }
    auto operator-(iterator_impl const& that) { return curr_pos - that.curr_pos; }
    auto& operator()() { return deref(particles); }
    auto& operator()() const { return deref(particles); }
    auto& operator*() _PHARE_ALL_FN_ { return *this; }
    auto& operator*() const _PHARE_ALL_FN_ { return *this; }

    auto idx() const { return curr_pos; }

    auto& weight() _PHARE_ALL_FN_ { return deref(particles).weight_[curr_pos]; }
    auto& weight() const _PHARE_ALL_FN_ { return deref(particles).weight_[curr_pos]; }
    auto& charge() _PHARE_ALL_FN_ { return deref(particles).charge_[curr_pos]; }
    auto& charge() const _PHARE_ALL_FN_ { return deref(particles).charge_[curr_pos]; }
    auto& iCell() _PHARE_ALL_FN_ { return deref(particles).iCell_[curr_pos]; }
    auto& iCell() const _PHARE_ALL_FN_ { return deref(particles).iCell_[curr_pos]; }
    auto& delta() _PHARE_ALL_FN_ { return deref(particles).delta_[curr_pos]; }
    auto& delta() const _PHARE_ALL_FN_ { return deref(particles).delta_[curr_pos]; }
    auto& v() _PHARE_ALL_FN_ { return deref(particles).v_[curr_pos]; }
    auto& v() const _PHARE_ALL_FN_ { return deref(particles).v_[curr_pos]; }

    Particle<dimension> copy() const _PHARE_ALL_FN_
    {
        return {weight(), charge(), iCell(), delta(), v()};
    }

    T particles;
    std::size_t curr_pos = 0;
};


template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using ParticleArray_SOAView = SoAParticles<SoASpan<dim, alloc_mode>>;


} // namespace PHARE::core


// namespace std
// {
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
//         std::copy(&src[src_begin.curr_pos], &src[src_end.curr_pos], &dst[dst_begin.curr_pos]);
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
//          typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template
//          iterator_impl<T>>
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
//          typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template
//          iterator_impl<T>>
// void iter_swap(Iter a, Iter b)
// {
//     PHARE_LOG_LINE_STR("");
//     PHARE_ASSERT(&a.particles == &b.particles);

//     a.particles.swap(a.curr_pos, b.curr_pos);
// }

// template<typename OuterSuper, typename T, typename Predicate,
//          typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template
//          iterator_impl<T>>
// Iter find_if_not(Iter first, Iter last, Predicate q)
// { // https://en.cppreference.com/w/cpp/algorithm/find
//     PHARE_LOG_LINE_STR("");
//     for (; first != last; ++first)
//         if (!q(first)) // pass iterator!
//             return first;
//     return last;
// }

// template<typename OuterSuper, typename T, typename Predicate/*,
//          typename Iter = typename PHARE::core::SoAParticles<OuterSuper>::template
//          iterator_impl<T>*/>
// typename PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T> partition(typename
// PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T> first, typename
// PHARE::core::SoAParticles<OuterSuper>::template iterator_impl<T> last, Predicate p)
// { // https://en.cppreference.com/w/cpp/algorithm/partition
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
// typename std::enable_if_t<Iterator<O>::is_soa, std::size_t> distance(Iterator<O> const& a,
//                                                                      Iterator<O> const& b)
// {
//     PHARE_ASSERT(&deref(a.particles).weight(0) == &deref(b.particles).weight(0));
//     PHARE_ASSERT(a.curr_pos <= b.curr_pos);

//     return b.curr_pos - a.curr_pos;
// }

// } // namespace std


namespace PHARE::core
{


template<typename Super>
SoAParticles<Super>::SoAParticles(SoAParticles<Super>::iterator start,
                                  SoAParticles<Super>::iterator end)
    : Super{std::distance(start, end)}
{
    std::copy(start, end, this->begin()); // impl above
}

} // namespace PHARE::core




#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP */
