#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SOA_HPP

#include "core/data/particles/particle.hpp"


#include "core/utilities/span.hpp" // for flatten

namespace PHARE::core
{
template<std::size_t dim, std::size_t size_>
struct SoAArray
{
    static constexpr auto dimension = dim;

    std::array<double, size_> weight_, charge_;
    std::array<std::array<int, dim>, size_> iCell_;
    std::array<std::array<double, dim>, size_> delta_;
    std::array<std::array<double, 3>, size_> v_;


    std::array<std::array<double, 3>, size_> E_;
    std::array<std::array<double, 3>, size_> B_;


    auto constexpr static size() { return size_; }
};




// used in SOA_ParticelArray::operator[](std::size_t) for std algorithm interactions
//  allows copying between vector and array impls without iterators
//   could be a better way maybe but it's not so simple to change
template<std::size_t dim>
using SoAParticle_crt = std::tuple<double const&,                  //  weight
                                   double const&,                  // charge
                                   std::array<int, dim> const&,    // iCell
                                   std::array<double, dim> const&, // delta
                                   std::array<double, 3> const&,   // v
                                   std::array<double, 3> const&,   // E
                                   std::array<double, 3> const&    // B
                                   >;



template<std::size_t dim>
struct SoAVector
{
    static constexpr auto dimension = dim;
    using container_type            = std::vector<Particle<dim>>;
    using value_type                = typename container_type::value_type;


    SoAVector() {}

    SoAVector(std::size_t size)
        : weight_(size)
        , charge_(size)
        , iCell_(size)
        , delta_(size)
        , v_(size)
        , E_(size)
        , B_(size)
    {
    }

    template<typename Particle_t>
    SoAVector(std::size_t size, Particle_t&& from)
        : weight_(size, from.weight())
        , charge_(size, from.charge())
        , iCell_(size, from.iCell())
        , delta_(size, from.delta())
        , v_(size, from.v())
        , E_(size, from.E())
        , B_(size, from.B())
    {
    }

    auto size() const { return weight_.size(); }

    std::vector<double> weight_, charge_;
    std::vector<std::array<int, dim>> iCell_;
    std::vector<std::array<double, dim>> delta_;
    std::vector<std::array<double, 3>> v_;

    std::vector<std::array<double, 3>> E_;
    std::vector<std::array<double, 3>> B_;
};


// used when the memory is owned elsewhere, e.g. numpy arrays
template<std::size_t dim, bool _const_ = false>
struct SoASpan
{
    static constexpr auto dimension = dim;


    template<typename Container_int, typename Container_double, bool const__ = _const_,
             typename = std::enable_if_t<const__>>
    SoASpan(Container_int const&& _iCell, Container_double const&& _delta,
            Container_double const&& _weight, Container_double const&& _charge,
            Container_double const&& _v)
        : size_{_weight.size()}
        , weight_{reinterpret_cast<double const*>(_weight.data())}
        , charge_{reinterpret_cast<double const*>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim> const*>(_iCell.data())}
        , delta_{reinterpret_cast<std::array<double, dim> const*>(_delta.data())}
        , v_{reinterpret_cast<std::array<double, 3> const*>(_v.data())}
    {
    }


    template<typename Container_int, typename Container_double, bool const__ = _const_,
             typename = std::enable_if_t<!const__>>
    SoASpan(Container_int&& _iCell, Container_double&& _delta, Container_double&& _weight,
            Container_double&& _charge, Container_double&& _v)
        : size_{_weight.size()}
        , weight_{reinterpret_cast<double* const>(_weight.data())}
        , charge_{reinterpret_cast<double* const>(_charge.data())}
        , iCell_{reinterpret_cast<std::array<int, dim>* const>(_iCell.data())}
        , delta_{reinterpret_cast<std::array<double, dim>* const>(_delta.data())}
        , v_{reinterpret_cast<std::array<double, 3>* const>(_v.data())}
    {
    }

    auto& size() const { return size_; }

    std::size_t size_;

    template<typename T>
    using container_t = std::conditional_t<_const_, T const*, T* const>;

    container_t<double> weight_, charge_;
    container_t<std::array<int, dim>> iCell_;
    container_t<std::array<double, dim>> delta_;
    container_t<std::array<double, 3>> v_;

    container_t<std::array<double, dim>> E_{};
    container_t<std::array<double, 3>> B_{};
};



template<typename Super_>
class SoAParticles : public Super_
{
    template<typename T>
    struct iterator_impl;

public:
    static constexpr bool is_contiguous = true;
    static constexpr auto dimension     = Super_::dimension;

    using Super = Super_;
    using This  = SoAParticles<Super>;
    using Super::size;
    using Particle_t = SoAParticle_crt<dimension>;
    using value_type = Particle_t;
    template<std::size_t size>
    using array_type = SoAParticles<SoAArray<dimension, size>>;

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    // public for pybind but avoid otherwise
    using Super::B_;
    using Super::charge_;
    using Super::delta_;
    using Super::E_;
    using Super::iCell_;
    using Super::v_;
    using Super::weight_;

protected:
    template<typename T>
    auto static constexpr is_vector()
    {
        return std::is_same_v<T, SoAVector<dimension>>;
    }
    template<typename T>
    auto static constexpr is_span()
    {
        return std::is_same_v<
                   T, SoASpan<dimension, true>> or std::is_same_v<T, SoASpan<dimension, false>>;
    }


public:
    SoAParticles() {}

    template<typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    SoAParticles(std::size_t size = 0)
        : Super{size}
    {
    }

    template<typename Particle_t, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    SoAParticles(std::size_t size, Particle_t&& particle)
        : Super{size, std::forward<Particle_t>(particle)}
    {
    }

    template<typename... Args, typename S = Super, typename = std::enable_if_t<is_span<S>()>>
    SoAParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }
    template<typename... Args, typename S = Super, typename = std::enable_if_t<is_span<S>()>>
    SoAParticles(Args const&... args)
        : Super{args...}
    {
    }

    SoAParticles(iterator start, iterator end);

    auto constexpr static size_of_particle()
    {
        return sizeof(typename decltype(iCell_)::value_type)
               + sizeof(typename decltype(delta_)::value_type)
               + sizeof(typename decltype(weight_)::value_type)
               + sizeof(typename decltype(charge_)::value_type)
               + sizeof(typename decltype(v_)::value_type)
               + sizeof(typename decltype(E_)::value_type)
               + sizeof(typename decltype(B_)::value_type);
    }



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


    auto& E(std::size_t i) const { return E_[i]; }
    auto& E(std::size_t i) { return E_[i]; }

    auto& B(std::size_t i) const { return B_[i]; }
    auto& B(std::size_t i) { return B_[i]; }

    auto& weight() const { return weight_; }
    auto& charge() const { return charge_; }
    auto& iCell() const { return iCell_; }
    auto& delta() const { return delta_; }
    auto& v() const { return v_; }


    // for performing the same operation across all vectors e.g. with std apply
    auto as_tuple(std::size_t i)
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i], this->E_[i], this->B_[i]);
    }
    auto as_tuple(std::size_t i) const
    {
        return std::forward_as_tuple(this->weight_[i], this->charge_[i], this->iCell_[i],
                                     this->delta_[i], this->v_[i], this->E_[i], this->B_[i]);
    }


    // to be avoided, but is convenient
    auto operator[](std::size_t i) { return as_tuple(i); }
    auto operator[](std::size_t i) const { return as_tuple(i); }



    auto begin() { return iterator(*this); }
    auto begin() const { return const_iterator(*this); }
    auto end() { return iterator(*this) + size(); }
    auto end() const { return const_iterator(*this) + size(); }



    template<typename Impl, typename S = Super, std::enable_if_t<is_vector<S>(), int> = 0>
    void push_back(iterator_impl<Impl> const& particle)
    {
        this->weight_.push_back(particle.weight());
        this->charge_.push_back(particle.charge());
        this->iCell_.push_back(particle.iCell());
        this->delta_.push_back(particle.delta());
        this->v_.push_back(particle.v());
        this->E_.push_back(particle.E());
        this->B_.push_back(particle.B());
    }


    template<typename Particle_t, typename S = Super, std::enable_if_t<is_vector<S>(), int> = 0>
    void push_back(Particle_t const& particle)
    {
        auto const& [w, c, i, d, v, E, B] = particle;
        this->weight_.push_back(w);
        this->charge_.push_back(c);
        this->iCell_.push_back(i);
        this->delta_.push_back(d);
        this->v_.push_back(v);

        this->E_.push_back(E);
        this->B_.push_back(B);
    }



    template<typename... Args>
    auto emplace_back(Args const&... args)
    {
        auto arg_tuple  = std::forward_as_tuple(args...);
        auto this_tuple = as_tuple<false>();
        for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto ic) {
            auto constexpr i = ic();
            std::get<i>(this_tuple).emplace_back(std::get<i>(arg_tuple));
        });
        this->E_.push_back(ConstArray<double, 3>());
        this->B_.push_back(ConstArray<double, 3>());

        return (*this)[this->size() - 1];
    }



    template<bool full = true>
    // full addes E/B to the tuple
    auto as_tuple()
    {
        if constexpr (full)
            return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_, E_, B_);
        else
            return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }

    template<bool full = true>
    // full addes E/B to the tuple
    auto as_tuple() const
    {
        if constexpr (full)
            return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_, E_, B_);
        else
            return std::forward_as_tuple(weight_, charge_, iCell_, delta_, v_);
    }


    void clear()
    {
        std::apply([](auto&... container) { ((container.clear()), ...); }, as_tuple());
    }

    void resize(std::size_t size)
    {
        std::apply([&](auto&... container) { ((container.resize(size)), ...); }, as_tuple());
    }

    void reserve(std::size_t size)
    {
        std::apply([&](auto&... container) { ((container.reserve(size)), ...); }, as_tuple());
    }

    void swap(This& that)
    {
        auto this_tuple = as_tuple();
        auto that_tuple = that.as_tuple();
        for_N<std::tuple_size_v<decltype(this_tuple)>>([&](auto ic) {
            auto constexpr i = ic();
            std::get<i>(this_tuple).swap(std::get<i>(that_tuple));
        });
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
        std::swap(E_[a], E_[b]);
        std::swap(B_[a], B_[b]);
    }



    void erase(iterator first, iterator last)
    {
        std::apply(
            [&](auto&... container) {
                ((container.erase(container.begin() + first.curr_pos,
                                  container.begin() + last.curr_pos)),
                 ...);
            },
            as_tuple());
    }

    struct Sorter;
    void sort() { Sorter{*this}(); }
    void sort(std::int64_t s, std::int64_t e) { Sorter{*this}(s, e); }
};


template<typename OuterSuper>
template<typename T>
struct SoAParticles<OuterSuper>::iterator_impl
{
    static constexpr auto dimension     = OuterSuper::dimension;
    auto static constexpr is_const      = std::is_const_v<T>;
    auto static constexpr is_contiguous = true;

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using value_type        = std::decay_t<T>;
    using pointer           = std::decay_t<T>*;
    using reference         = std::decay_t<T>&;


    iterator_impl(T& particles_)
        : particles{particles_}
    {
    }

    iterator_impl(iterator_impl&& that)      = default;
    iterator_impl(iterator_impl const& that) = default;

    auto& operator++()
    {
        ++curr_pos;
        return *this;
    }
    auto& operator+=(std::int64_t i)
    {
        curr_pos += i;
        return *this;
    }

    auto& operator--()
    {
        --curr_pos;
        return *this;
    }
    auto operator+(std::int64_t i)
    {
        auto copy = *this;
        copy.curr_pos += i;
        return copy;
    }
    auto operator-(std::int64_t i)
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

    auto& operator()() { return particles; }
    auto& operator()() const { return particles; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

    auto idx() const { return curr_pos; }


    auto& weight() { return particles.weight_[curr_pos]; }
    auto& weight() const { return particles.weight_[curr_pos]; }
    auto& charge() { return particles.charge_[curr_pos]; }
    auto& charge() const { return particles.charge_[curr_pos]; }
    auto& iCell() { return particles.iCell_[curr_pos]; }
    auto& iCell() const { return particles.iCell_[curr_pos]; }
    auto& delta() { return particles.delta_[curr_pos]; }
    auto& delta() const { return particles.delta_[curr_pos]; }
    auto& v() { return particles.v_[curr_pos]; }
    auto& v() const { return particles.v_[curr_pos]; }

    auto& E() { return particles.E_[curr_pos]; }
    auto& E() const { return particles.E_[curr_pos]; }
    auto& B() { return particles.B_[curr_pos]; }
    auto& B() const { return particles.B_[curr_pos]; }

    T& particles;
    std::size_t curr_pos = 0;
};



template<typename Super>
struct SoAParticles<Super>::Sorter
{
    template<typename V>
    auto el_wise_less(V const& v0, V const& v1)
    {
        for (std::int16_t i = 0; i < v0.size(); ++i)
            if (v0[i] < v1[i])
                return true;
            else if (v0[i] != v1[i])
                return false;
        return false;
    }
    template<typename V>
    auto el_wise_gr8r(V const& v0, V const& v1)
    {
        for (std::int16_t i = 0; i < v0.size(); ++i)
            if (v0[i] > v1[i])
                return true;
            else if (v0[i] != v1[i])
                return false;
        return false;
    }

    void operator()(std::int64_t l, std::int64_t r)
    {
        auto i = l;
        auto j = r;

        auto const& iCells = self.iCell_;
        auto const half    = iCells[(l + r) / 2];

        do
        {
            while (el_wise_less(iCells[i], half))
                i++;
            while (el_wise_gr8r(iCells[j], half))
                j--;

            if (i <= j)
            {
                self.swap(i, j);

                i++;
                j--;
            }
        } while (i <= j);

        if (l < j)
            (*this)(l, j);
        if (i < r)
            (*this)(i, r);
    }

    void operator()() { (*this)(0, self.size() - 1); }

    SoAParticles<Super>& self;
};


template<std::size_t dim, bool _const_ = false>
using ParticleArray_SOAView = SoAParticles<SoASpan<dim, _const_>>;


} // namespace PHARE::core


namespace std
{
using namespace PHARE::core;


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, Particle<Iterator<O>::dimension>>
copy(Iterator<O> src)
{
    return {src.weight(), src.charge(), src.iCell(), src.delta(), src.v()};
}


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void> copy(Iterator<O> src_begin, //
                                                                 Iterator<O> src_end,
                                                                 Iterator<O> dst_begin)
{
    auto src_tuple = src_begin.particles.as_tuple();
    auto dst_tuple = dst_begin.particles.as_tuple();
    for_N<std::tuple_size_v<decltype(src_tuple)>>([&](auto ic) {
        auto constexpr i = ic();
        auto& src        = std::get<i>(src_tuple);
        auto& dst        = std::get<i>(dst_tuple);
        std::copy(&src[src_begin.curr_pos], &src[src_end.curr_pos], &dst[dst_begin.curr_pos]);
    });
}



template<template<typename> typename Iterator, typename O, typename Inserter>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void> copy(Iterator<O> src_begin, //
                                                                 Iterator<O> src_end, Inserter i)
{
    auto const& particles = src_begin.particles;

    for (; src_begin != src_end; ++src_begin, ++i)
        *i = particles[src_begin.curr_pos];
}

template<template<typename> typename Iterator, typename O, typename Inserter, typename Fn>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void>
copy_if(Iterator<O> src_begin, Iterator<O> src_end, Inserter i, Fn check)
{
    auto const& particles = src_begin.particles;

    for (; src_begin != src_end; ++src_begin)
        if (check(src_begin))
        {
            *i = particles[src_begin.curr_pos];
            ++i;
        }
}


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void> iter_swap(Iterator<O>& a,
                                                                      Iterator<O>& b)
{
    assert(&a.particles == &b.particles);

    a.particles.swap(a.curr_pos, b.curr_pos);
}

template<template<typename> typename Iterator, typename O, typename Predicate>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, Iterator<O>>
find_if_not(Iterator<O> first, Iterator<O> last, Predicate q)
{ // https://en.cppreference.com/w/cpp/algorithm/find
    for (; first != last; ++first)
        if (!q(first)) // pass iterator!
            return first;
    return last;
}


template<template<typename> typename Iterator, typename O, typename Predicate>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, Iterator<O>>
partition(Iterator<O> first, Iterator<O> last, Predicate p)
{ // https://en.cppreference.com/w/cpp/algorithm/partition
    first = std::find_if_not(first, last, p);
    if (first == last)
        return first;

    for (auto i = first + 1; i != last; ++i)
    {
        if (p(i)) // pass iterator!
        {
            std::iter_swap(i, first);
            ++first;
        }
    }
    return first;
}


template<template<typename> typename Iterator, typename O, typename Fn>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, void>
transform(Iterator<O> src_begin, Iterator<O> src_end, Iterator<O> dst_begin, Fn transform)
{
    if (&src_begin.particles != &dst_begin.particles)
        throw std::runtime_error("transform cannot work");

    for (; src_begin != src_end; ++src_begin, ++dst_begin)
        /*dst_begin = */ transform(src_begin); // transform edits in place
}


template<template<typename> typename Iterator, typename O>
// SFINAE to support const iterators
typename std::enable_if_t<Iterator<O>::is_contiguous, std::size_t> distance(Iterator<O> const& a,
                                                                            Iterator<O> const& b)
{
    assert(&a.particles == &b.particles);
    assert(a.curr_pos <= b.curr_pos);

    return b.curr_pos - a.curr_pos;
}



template<typename Super>
void sort(SoAParticles<Super>& particles)
{
    particles.sort();
}


} // namespace std



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
