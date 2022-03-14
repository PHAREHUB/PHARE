#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP


#include <vector>
#include <cstddef>
#include <utility>

#include "particle.hpp"

#include "core/utilities/iterators.hpp"
#include "core/utilities/point/point.hpp"

namespace PHARE::core
{
template<std::size_t dim, std::size_t size_>
class AoSArray
{
public:
    static constexpr auto dimension = dim;
    using Particle_t                = Particle<dim>;
    using container_type            = std::array<Particle_t, size_>;


    using iterator       = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

    auto begin() const { return particles_.begin(); }
    auto begin() { return particles_.begin(); }

    auto end() const { return particles_.end(); }
    auto end() { return particles_.end(); }

    auto constexpr static size() { return size_; }

protected:
    container_type particles_;
};

template<std::size_t dim>
class AoSVector
{
    template<typename Iterator>
    auto check_distance_size_t(Iterator const& start, Iterator const& end)
    {
        auto dist = std::distance(start, end);
        if (dist < 0)
            throw std::runtime_error("Error, number must be postive");
        return static_cast<std::size_t>(dist);
    }

    using This = AoSVector<dim>;

public:
    static constexpr auto dimension = dim;

    using Particle_t     = Particle<dim>;
    using container_type = std::vector<Particle_t>;
    using value_type     = Particle_t;

    AoSVector(std::size_t size = 0)
        : particles_(size)
    {
    }

    template<typename Particle_t>
    AoSVector(std::size_t size, Particle_t&& particle)
        : particles_(size, particle)
    {
    }

    template<typename Iterator>
    AoSVector(Iterator start, Iterator end)
        : AoSVector{check_distance_size_t(start, end)}
    {
        std::copy(start, end, particles_.begin());
    }

    template<typename T>
    struct iterator_impl;

    using iterator       = iterator_impl<This>;
    using const_iterator = iterator_impl<This const>;

    auto size() const { return particles_.size(); }

    void clear() { return particles_.clear(); }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }

    auto& operator[](std::size_t i) const { return particles_.data()[i]; }
    auto& operator[](std::size_t i) { return particles_.data()[i]; }

    bool operator==(This const& that) const { return (this->particles_ == that.particles_); }

    auto begin() const { return const_iterator{particles_.begin(), *this}; }
    auto begin() { return iterator{particles_.begin(), *this}; }

    auto end() const { return const_iterator{particles_.end(), *this}; }
    auto end() { return iterator{particles_.end(), *this}; }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles_.insert(position, first, last);
    }

    auto back() { return particles_.back(); }
    auto front() { return particles_.front(); }

    iterator erase(iterator position)
    {
        return {particles_.erase(position), *this};
        ;
    }

    template<typename It0, typename It1>
    iterator erase(It0 first, It1 last)
    {
        return {particles_.erase(first, last), *this};
    }

    Particle_t& emplace_back() { return particles_.emplace_back(); }
    Particle_t& emplace_back(Particle_t&& p) { return particles_.emplace_back(p); }
    Particle_t& emplace_back(Particle_t const& p) { return particles_.emplace_back(p); }

    template<typename... Args>
    Particle_t& emplace_back(Args const&... args)
    {
        return particles_.emplace_back(args...);
    }
    template<typename... Args>
    Particle_t& emplace_back(Args&&... args)
    {
        return particles_.emplace_back(Particle_t{args...});
    }

    void push_back(Particle_t const& p) { particles_.push_back(p); }
    void push_back(Particle_t&& p) { particles_.push_back(std::forward<Particle_t>(p)); }

    template<typename That>
    void swap(That& that)
    {
        std::swap(this->particles_, that.particles_);
    }



protected:
    container_type particles_;
};

template<bool is_vector>
struct AoSParticleIterator
{
};

template<typename Super_>
struct AoSParticles : public Super_
{
    using Super = Super_;
    using This  = AoSParticles<Super>;

    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = Super::dimension;

    using Super::particles_;
    using container_type = typename Super::container_type;
    using Particle_t     = typename Super::Particle_t;

    template<std::size_t size>
    using array_type = AoSParticles<AoSArray<dimension, size>>;

    template<typename T>
    auto static constexpr is_vector()
    {
        return std::is_same_v<T, AoSVector<dimension>>;
    }

    AoSParticles() {}

    template<typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(std::size_t size)
        : Super(size)
    {
    }

    template<typename Particle_t, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(std::size_t size, Particle_t&& particle)
        : Super(size, std::forward<Particle_t>(particle))
    {
    }


    template<typename Iterator, typename S = Super, typename = std::enable_if_t<is_vector<S>()>>
    AoSParticles(Iterator start, Iterator end)
        : Super{start, end}
    {
    }


    template<bool is_const = false> // :(
    auto static constexpr iterator_type()
    {
        throw std::runtime_error("never to be called in non constexpr context");
        if constexpr (is_vector<Super>())
            return std::conditional_t<is_const, typename Super::template iterator_impl<This const>*,
                                      typename Super::template iterator_impl<This>*>{nullptr};
        else
            return std::conditional_t<is_const, typename Super::const_iterator,
                                      typename Super::iterator>{nullptr};
    }

    using iterator       = std::decay_t<decltype(*iterator_type())>;
    using const_iterator = std::decay_t<decltype(*iterator_type<true>())>;

    auto begin() const
    {
        if constexpr (is_vector<Super>())
            return const_iterator{particles_.begin(), *this};
        else
            return Super::begin();
    }
    auto begin()
    {
        if constexpr (is_vector<Super>())
            return iterator{particles_.begin(), *this};
        else
            return Super::begin();
    }

    auto end() const
    {
        if constexpr (is_vector<Super>())
            return const_iterator{particles_.end(), *this};
        else
            return Super::end();
    }
    auto end()
    {
        if constexpr (is_vector<Super>())
            return iterator{particles_.end(), *this};
        else
            return Super::end();
    }



    auto& weight(std::size_t i) const { return particles_[i].weight(); }
    auto& weight(std::size_t i) { return particles_[i].weight(); }

    template<typename Weight>
    void weight(std::size_t i, Weight weight)
    {
        particles_[i].weight() = weight;
    }

    auto& charge(std::size_t i) const { return particles_[i].charge(); }
    auto& charge(std::size_t i) { return particles_[i].charge(); }

    template<typename Charge>
    void charge(std::size_t i, Charge charge)
    {
        particles_[i].charge() = charge;
    }

    auto& iCell(std::size_t i) const { return particles_[i].iCell(); }
    auto& iCell(std::size_t i) { return particles_[i].iCell(); }

    template<typename ICell>
    void iCell(std::size_t i, ICell iCell)
    {
        particles_[i].iCell() = iCell;
    }

    auto& delta(std::size_t i) const { return particles_[i].delta(); }
    auto& delta(std::size_t i) { return particles_[i].delta(); }

    template<typename Delta>
    void delta(std::size_t i, Delta delta)
    {
        particles_[i].delta() = delta;
    }

    auto& v(std::size_t i) const { return particles_[i].v(); }
    auto& v(std::size_t i) { return particles_[i].v(); }

    template<typename V>
    void v(std::size_t i, V v)
    {
        particles_[i].v() = v;
    }

    std::size_t size() const { return particles_.size(); }

    template<typename S = Super, std::enable_if_t<is_vector<S>(), bool> = 0>
    std::size_t capacity() const
    {
        return particles_.capacity();
    }


    auto& E(std::size_t i) const { return particles_[i].E_; }
    auto& E(std::size_t i) { return particles_[i].E_; }

    auto& B(std::size_t i) const { return particles_[i].B_; }
    auto& B(std::size_t i) { return particles_[i].B_; }


    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }

    auto& vector() { return particles_; }
    auto& vector() const { return particles_; }

    auto constexpr static size_of_particle() { return sizeof(Particle_t); }
};



template<std::size_t dim>
template<typename T>
struct AoSVector<dim>::iterator_impl
    : public wrapped_iterator<T, typename std::decay_t<T>::container_type>
{
    static constexpr auto is_contiguous = false;
    static constexpr auto dimension     = dim;

    using Super = wrapped_iterator<T, typename std::decay_t<T>::container_type>;

    iterator_impl() {}

    template<typename Iterator>
    iterator_impl(Iterator iter, T& container)
        : Super{iter, &container}
    {
    }
    iterator_impl(iterator_impl const& that) = default;

    iterator_impl& operator=(iterator_impl const& other) = default;

    bool operator==(iterator_impl const& other) const
    {
        bool superbool = (static_cast<Super const&>(*this) == static_cast<Super const&>(other));
        return superbool;
    }

    auto operator+(std::size_t i)
    {
        iterator_impl copy = *this;
        static_cast<Super&>(copy) += i;
        return copy;
    }

    auto& weight() { return (*this)->weight(); }
    auto& weight() const { return (*this)->weight(); }
    auto& charge() { return (*this)->charge(); }
    auto& charge() const { return (*this)->charge(); }
    auto& iCell() { return (*this)->iCell(); }
    auto& iCell() const { return (*this)->iCell(); }
    auto& delta() { return (*this)->delta(); }
    auto& delta() const { return (*this)->delta(); }
    auto& v() { return (*this)->v(); }
    auto& v() const { return (*this)->v(); }

    template<typename Weight>
    void weight(Weight weight)
    {
        (*this)->weight_ = weight;
    }

    template<typename Charge>
    void charge(Charge charge)
    {
        (*this)->charge_ = charge;
    }

    template<typename ICell>
    void iCell(ICell iCell)
    {
        (*this)->iCell_ = iCell;
    }

    template<typename Delta>
    void delta(Delta delta)
    {
        (*this)->delta_ = delta;
    }

    template<typename V>
    void v(V v)
    {
        (*this)->v_ = v;
    }
};

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP */
