#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP

#include "core/data/vector.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

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
    void assign(Particles_t const& src, std::size_t const idx, std::size_t const dst)
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
    void assign(Particle_t const& src, std::size_t const dst)
    {
        assert(dst < size_);
        this->weight_[dst] = src.weight();
        this->charge_[dst] = src.charge();
        this->iCell_[dst]  = src.iCell();
        this->delta_[dst]  = src.delta();
        this->v_[dst]      = src.v();
    }

    auto as_tuple() { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
    auto as_tuple() const { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
};

// used when the memory is owned elsewhere, e.g. numpy arrays
template<std::size_t dim, auto alloc_mode_, auto _const_ = 0>
struct SoASpan
{
    static_assert(std::is_same_v<decltype(alloc_mode_), AllocatorMode>);

    auto static constexpr alloc_mode   = alloc_mode_;
    auto static constexpr storage_mode = StorageMode::SPAN;
    auto static constexpr dimension    = dim;
    template<typename T>
    using container_t = std::conditional_t<_const_, T const*, T*>;

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

    auto size() const { return size_; }
    void clear() { size_ = 0; }
    void resize(std::size_t const& s)
    {
        PHARE_ASSERT(s <= size_); // can't be bigger
        size_ = s;
    }

    void pop_back() { --size_; }
    auto size_address() { return &size_; }

    auto as_tuple() { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
    auto as_tuple() const { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }

    template<typename Particles_t>
    void reset(Particles_t& particles)
    {
        auto self = as_tuple();
        auto that = particles.as_tuple();
        for_N<std::tuple_size_v<decltype(that)>>(
            [&](auto vi) { std::get<vi>(self) = std::get<vi>(that).data(); });
        size_ = particles.size();
    }

    std::size_t size_;
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

    template<typename Type>
    using container_t = std::vector<Type>;

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

    auto size() const { return weight_.size(); }

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
        std::apply([&](auto&... container) { ((container.resize(size)), ...); }, as_tuple());
    }

    template<auto type>
    void on_moved()
    {
        // noop
    }

    template<auto type>
    void on_appended()
    {
        // noop
    }

    container_t<double> weight_, charge_;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta_;
    container_t<std::array<double, 3>> v_;

    template<typename V>
    static auto& get_vec(V& v) // we can probably delete this now
    {
        return v;
    }
};


template<typename Super_>
class SoAParticles : public Super_
{
public:
    using Super                        = Super_;
    using This                         = SoAParticles<Super>;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr layout_mode  = LayoutMode::SoA;
    auto static constexpr storage_mode = Super::storage_mode;

    using Particle_t = SoAParticle_crt<dimension>;
    using Super::size;

    using Span_t = SoAParticles<SoASpan<dimension, alloc_mode>>;
    friend class SoAParticles<SoASpan<dimension, alloc_mode>>;

    // public for pybind but avoid otherwise
    using Super::charge_;
    using Super::delta_;
    using Super::iCell_;
    using Super::v_;
    using Super::weight_;

    template<typename... Args>
    SoAParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
        : Super{std::forward<Args>(args)...}
    {
    }

    SoAParticles(This const& that)    = default;
    SoAParticles(This&& that)         = default;
    This& operator=(This&& that)      = default;
    This& operator=(This const& that) = default;

    // no begin()/end() here - SoA has no single Particle_t& to hand back (fields live in
    // separate arrays), only index-based access (see weight(i)/charge(i)/... below).

    auto& weight(std::size_t i) const { return weight_[i]; }
    auto& weight(std::size_t i) { return weight_[i]; }

    auto& charge(std::size_t i) const { return charge_[i]; }
    auto& charge(std::size_t i) { return charge_[i]; }

    auto& iCell(std::size_t i) const { return iCell_[i]; }
    auto& iCell(std::size_t i) { return iCell_[i]; }

    auto& delta(std::size_t i) const { return delta_[i]; }
    auto& delta(std::size_t i) { return delta_[i]; }

    auto& v(std::size_t i) const { return v_[i]; }
    auto& v(std::size_t i) { return v_[i]; }

    auto& weight() { return weight_; }
    auto& charge() { return charge_; }
    auto& iCell() { return iCell_; }
    auto& delta() { return delta_; }
    auto& v() { return v_; }

    auto& weight() const { return weight_; }
    auto& charge() const { return charge_; }
    auto& iCell() const { return iCell_; }
    auto& delta() const { return delta_; }
    auto& v() const { return v_; }

    // for performing the same operation across all vectors e.g. with std apply
    auto as_tuple(std::size_t i)
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i]);
    }

    auto as_tuple(std::size_t i) const
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i]);
    }

    auto as_tuple() { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }
    auto as_tuple() const { return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_); }

    template<typename Particle_t, auto S = storage_mode>
        requires(S == StorageMode::VECTOR)
    void push_back(Particle_t const& particle)
    {
        auto const& [w, c, i, d, v] = particle;
        Super::get_vec(this->weight_).push_back(w);
        Super::get_vec(this->charge_).push_back(c);
        Super::get_vec(this->iCell_).push_back(i);
        Super::get_vec(this->delta_).push_back(d);
        Super::get_vec(this->v_).push_back(v);
    }


    template<typename Particle_t, auto S = storage_mode>
        requires(S == StorageMode::VECTOR)
    auto emplace_back(Particle_t const& particle)
    {
        Super::get_vec(this->weight_).emplace_back(particle.weight());
        Super::get_vec(this->charge_).emplace_back(particle.charge());
        Super::get_vec(this->iCell_).emplace_back(particle.iCell());
        Super::get_vec(this->delta_).emplace_back(particle.delta());
        Super::get_vec(this->v_).emplace_back(particle.v());
        return size() - 1; // index of the appended particle - SoA has no per-particle
                           // reference/iterator to hand back, unlike std::vector
    }

    template<typename That, auto S = storage_mode,
             typename
             = std::enable_if_t<S == StorageMode::VECTOR and That::layout_mode == LayoutMode::SoA>>
    void emplace_back(That const& src);


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
    auto emplace_back(Args const&... args);


    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    void reserve(std::size_t const& size)
    {
        std::apply([&](auto&... container) { ((container.reserve(size)), ...); }, as_tuple());
    }

    void swap(This& that)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = that.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto i) { std::get<i>(this_tuple).swap(std::get<i>(that_tuple)); });
    }

    void swap(std::size_t const& a, std::size_t const& b);

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    std::size_t capacity() const
    {
        return weight_.capacity(); // they're all the same
    }

    auto copy(std::size_t i) const
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

    void assign(std::size_t const& src, std::size_t const& dst)
    {
        std::apply([&](auto&... v) { ((v[dst] = v[src]), ...); }, as_tuple());
    }

    template<typename Particle_t>
    void assign(Particle_t const& src, std::size_t const& dst)
    {
        auto this_tuple = as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto i) { std::get<i>(this_tuple)[dst] = std::get<i>(*src); });
    }

    template<typename _Particles>
    void assign(_Particles const& src, std::size_t const& idx, std::size_t const& dst)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = src.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>(
            [&](auto vi) { std::get<vi>(this_tuple)[dst] = std::get<vi>(that_tuple)[idx]; });
    }
};

template<typename Super_>
template<typename That, auto S, typename>
void SoAParticles<Super_>::emplace_back(That const& src)
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

template<typename Super_>
template<typename... Args>
auto SoAParticles<Super_>::emplace_back(Args const&... args)
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

    return size() - 1; // index of the appended particle - see emplace_back(Particle_t const&)
}

template<typename Super_>
void SoAParticles<Super_>::swap(std::size_t const& a, std::size_t const& b)
{
    if (a == b)
        return;

    std::swap(weight_[a], weight_[b]);
    std::swap(charge_[a], charge_[b]);
    std::swap(iCell_[a], iCell_[b]);
    std::swap(delta_[a], delta_[b]);
    std::swap(v_[a], v_[b]);
}


// Per-index particle-like proxy for generic code (e.g. amr::Splitter) written against a
// Particle_t interface (.weight()/.charge()/.iCell()/.delta()/.v()) that SoA can't satisfy
// with a real begin()/end()/operator[] (see "no begin()/end() here" above). Reads and writes
// go straight through to the underlying SoA arrays at index `i`, no copy involved.
template<typename SoAParticles_t>
struct SoAParticleView
{
    static constexpr auto dimension = SoAParticles_t::dimension;

    auto& weight() const { return particles.weight(i); }
    auto& charge() const { return particles.charge(i); }
    auto& iCell() const { return particles.iCell(i); }
    auto& delta() const { return particles.delta(i); }
    auto& v() const { return particles.v(i); }

    SoAParticles_t& particles;
    std::size_t i;
};


template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using SoAVectorParticles = SoAParticles<SoAVector<dim, alloc_mode>>;

template<std::size_t dim, std::size_t size>
using SoAArrayParticles = SoAParticles<SoAArray<dim, size>>;

template<std::size_t dim, auto alloc_mode = AllocatorMode::CPU>
using ParticleArray_SOAView = SoAParticles<SoASpan<dim, alloc_mode>>;


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP */
