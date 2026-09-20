#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP

#include "core/def/phare_config.hpp" // IWYU pragma: keep

#include "core/data/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_type_options.hpp" // IWYU pragma: keep

#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace PHARE::core
{

template<auto opts>
class AoSArray
{
public:
    auto static constexpr dimension    = opts.dim;
    auto static constexpr alloc_mode   = AllocatorMode::CPU;
    auto static constexpr storage_mode = StorageMode::ARRAY;
    using Particle_t                   = ParticleDefaults<opts.dim>::Particle_t;
    using container_type               = std::array<Particle_t, opts.size>;

    auto begin() const { return particles_.begin(); }
    auto begin() { return particles_.begin(); }
    auto end() const { return particles_.end(); }
    auto end() { return particles_.end(); }
    auto constexpr static size() { return opts.size; }
    auto& operator[](std::size_t const i) const { return particles_.data()[i]; }
    auto& operator[](std::size_t const i) { return particles_.data()[i]; }

protected:
    container_type particles_;
};



template<auto opts>
class AoSSpan
{
    using This = AoSSpan<opts>;

public:
    auto static constexpr dimension    = opts.dim;
    auto static constexpr alloc_mode   = opts.alloc_mode;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = ParticleDefaults<opts.dim>::Particle_t;

    AoSSpan() = default;

    template<typename Container>
    AoSSpan(Container& container)
        : particles_{container.data(), container.size()}
    {
    }

    template<typename Container>
    AoSSpan(Container& container, std::size_t const& beg, std::size_t const& siz)
        : particles_{container.data() + beg, siz}
    {
    }
    template<typename Container>
    AoSSpan(Container& array, std::size_t const& siz)
        : AoSSpan{array, 0, siz}
    {
    }


    auto size() { return particles_.size(); }
    void clear() { particles_.s = 0; }
    void resize(std::size_t const& s) // potential out of bounds - so be careful!
    {
        particles_.s = s;
    }

    auto& operator[](std::size_t const i) const { return particles_[i]; }
    auto& operator[](std::size_t const i) { return particles_[i]; }


    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }

    void pop_back() { --particles_.s; }

    auto size_address() { return &particles_.s; }


    template<typename Particles_t>
    void reset(Particles_t& particles)
    {
        particles_.ptr = particles.data();
        particles_.s   = particles.size();
    }



    Span<Particle_t> particles_;
};


template<auto opts>
class AoSMappedSpan : public AoSSpan<opts>
{
    using box_t     = Box<int, opts.dim>;
    using CellMap_t = CellMap<CellMapOptions{opts.dim, opts.alloc_mode, StorageMode::SPAN}>;

public:
    using Super = AoSSpan<opts>;


    template<typename Particles_t>
    AoSMappedSpan(Particles_t& particles)
        : Super{particles}
        , box_{particles.box()}
        , cellMap_{particles.cellMap_}
    {
    }


    // auto map_size_address()  { return &particles_.s; }

    box_t box_;
    CellMap_t cellMap_;
};



template<auto opts>
class AoSVector
{
    using This = AoSVector<opts>;


    template<typename Iterator>
    auto check_distance_size_t(Iterator const& start, Iterator const& end)
    {
        auto dist = std::distance(start, end);
        if (dist < 0)
            throw std::runtime_error("Error, number must be postive");
        return static_cast<std::size_t>(dist);
    }

public:
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr alloc_mode   = opts.alloc_mode;
    auto static constexpr dimension    = opts.dim;

    using Particle_t     = ParticleDefaults<opts.dim>::Particle_t;
    using value_type     = Particle_t;
    using container_type = std::vector<Particle_t>;


    AoSVector(AoSVector&& that)
        : particles_(std::move(that.particles_))
    {
    }

    AoSVector(AoSVector const& from)            = default;
    AoSVector& operator=(AoSVector&& from)      = default;
    AoSVector& operator=(AoSVector const& from) = default;

    AoSVector(std::size_t size = 0)
        : particles_(size)
    {
    }

    template<typename Particle_t>
    AoSVector(std::size_t size, Particle_t const& particle)
        : particles_(size, particle)
    {
    }

    template<typename Iterator>
    AoSVector(Iterator start, Iterator end)
        : AoSVector{check_distance_size_t(start, end)}
    {
        std::copy(start, end, particles_.begin());
    }




    auto size() const { return particles_.size(); }
    auto capacity() const { return particles_.capacity(); }
    void clear() { particles_.clear(); }
    void reserve(std::size_t newSize) { particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }


    template<auto type>
    void on_moved() // noop
    {
    }

    template<auto type>
    void on_appended() // noop
    {
    }

    auto& operator[](std::size_t const i) const { return particles_.data()[i]; }
    auto& operator[](std::size_t const i) { return particles_.data()[i]; }
    bool operator==(This const& that) const { return (this->particles_ == that.particles_); }

    void pop_back() { particles_.pop_back(); }

    auto back() { return particles_.back(); }
    auto front() { return particles_.front(); }

    template<typename Iterator>
    auto erase(Iterator first, Iterator last);

    Particle_t& emplace_back() { return get_vec(particles_).emplace_back(); }

    Particle_t& emplace_back(Particle_t&& p)
    {
        return get_vec(particles_).emplace_back(std::forward<Particle_t>(p));
    }

    Particle_t& emplace_back(Particle_t const& p) { return get_vec(particles_).emplace_back(p); }

    template<typename Src>
    auto& emplace_back(Src const& src, std::size_t const& idx)
    {
        return emplace_back(src[idx]);
    }

    void emplace_back(This const& src)
    {
        for (auto const& p : src.particles_)
            emplace_back(p);
    }

    template<typename... Args>
    Particle_t& emplace_back(double const& weight, Args const&... args)
    {
        return get_vec(particles_).emplace_back(weight, args...);
    }
    template<typename... Args>
    Particle_t& emplace_back(double const& weight, Args&&... args)
    {
        return get_vec(particles_).emplace_back(weight, args...);
    }
    template<typename... Args>
    Particle_t& emplace_back(Args&&... args)
    {
        return get_vec(particles_).emplace_back(args...);
    }


    template<typename Src>
    void append(Src const& src, std::size_t const start, std::size_t const siz)
    {
        reserve(size() + siz);

        for (std::size_t i = 0; i < siz; ++i)
            emplace_back(src[i + start]);
    }

    void push_back(Particle_t const& p) { get_vec(particles_).push_back(p); }
    void push_back(Particle_t&& p) { get_vec(particles_).push_back(std::forward<Particle_t>(p)); }

    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }

protected:
    container_type particles_;

    template<typename V>
    static auto& get_vec(V& v)
    {
        return v;
    }
};

template<auto opts>
template<typename Iterator>
auto AoSVector<opts>::erase(Iterator first, Iterator last)
{
    // should we erase particles indexes associated with these iterators from the cellmap?
    // probably it does not matter if not. The reason is that
    // particles erased from the particlearray are so because they left
    // the patch cells to an outside cell.
    // But in principle that cell will never be accessed because it is outside the patch.
    // The only thing "bad" if these indexes are not deleted is that the
    // size of the cellmap becomes unequal to the size of the particleArray.
    // but  ¯\_(ツ)_/¯
    return particles_.erase(particles_.begin() + first.curr_pos,
                            particles_.begin() + last.curr_pos);
}


template<auto opts>
class AoSMappedVector : public AoSVector<opts>
{
    using box_t     = Box<int, opts.dim>;
    using CellMap_t = CellMap<CellMapOptions{opts.dim, opts.alloc_mode, StorageMode::VECTOR}>;

public:
    using Super = AoSVector<opts>;
    using Super::erase;

protected:
    using container_type = Super::container_type;
    using iterator_type  = container_type::iterator;

public:
    AoSMappedVector(box_t const& box = {}, auto&&... args)
        : Super{args...}
        , box_{box}
        , cellMap_{grow(box_, 1)} // safe outer shell == +1
    {
    }

    void erase(box_t const& box) { cellMap_.erase(this->particles_, box); }

    auto& map() { return cellMap_; }

protected:
    box_t box_;
    mutable CellMap_t cellMap_;
};


template<template<auto> typename Super_, auto opts>
struct AoSParticles : public Super_<opts>
{
    using Super      = Super_<opts>;
    using This       = AoSParticles;
    using Particle_t = typename Super::Particle_t;

    auto static constexpr dimension    = Super::dimension;
    auto static constexpr layout_mode  = LayoutMode::AoS;
    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Span_t = AoSParticles<AoSSpan, opts>;
    using Super::particles_;

    AoSParticles(AoSParticles&& that)
        : Super{std::forward<AoSParticles>(that)}
    {
    }
    AoSParticles(AoSParticles const&)            = default;
    AoSParticles& operator=(AoSParticles&&)      = default;
    AoSParticles& operator=(AoSParticles const&) = default;

    template<typename... Args>
    AoSParticles(Args&&... args)
        requires(self_excluding_constructible<This, Super, Args...>())
        : Super{std::forward<Args>(args)...}
    {
    }


    template<typename T>
    struct iterator_impl;

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto erase(iterator_impl<This*> a, iterator_impl<This*> b)
    {
        return Super::erase(a, b);
    }

    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args)
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T*>{t, args...};
    }
    auto begin() const { return it(this); }
    auto begin() { return it(this); }
    auto end() const { return it(this, size()); }
    auto end() { return it(this, size()); }


    auto& weight(std::size_t const i) const { return particles_[i].weight(); }
    auto& weight(std::size_t const i) { return particles_[i].weight(); }

    auto& charge(std::size_t const i) const { return particles_[i].charge(); }
    auto& charge(std::size_t const i) { return particles_[i].charge(); }

    auto& iCell(std::size_t const i) const { return particles_[i].iCell(); }
    auto& iCell(std::size_t const i) { return particles_[i].iCell(); }

    auto& delta(std::size_t const i) const { return particles_[i].delta(); }
    auto& delta(std::size_t const i) { return particles_[i].delta(); }

    auto& v(std::size_t const i) const { return particles_[i].v(); }
    auto& v(std::size_t const i) { return particles_[i].v(); }

    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto capacity() const
    {
        return particles_.capacity();
    }
    auto size() const { return particles_.size(); }

    template<typename IndexRange, typename Predicate>
    auto partition(IndexRange&& range, Predicate&& pred)
    {
        return std::partition(range.begin(), range.end(), pred);
    }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto& vector()
    {
        return particles_;
    }
    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto& vector() const
    {
        return particles_;
    }

    // does not swap cellmap (never did)
    void swap(This& that) { std::swap(this->particles_, that.particles_); }
    void swap(std::size_t const& a, std::size_t const& b)
    {
        if (a == b)
            return;
        std::swap(particles_[a], particles_[b]);
    }

    void check() const {}

    template<typename _Particles>
    void assign(_Particles const& src, std::size_t const& idx, std::size_t const& dst)
    {
        particles_[dst] = src[idx];
    }

    void assign(Particle_t const& src, std::size_t const& dst) { particles_[dst] = src; }

    void assign(std::size_t const& src, std::size_t const& dst)
    {
        particles_[dst] = particles_[src];
    }

    auto nbr_particles_in(Box<int, dimension> const /*box*/) const
    {
        throw std::runtime_error("finish this");

        return 10;
    }
};



template<template<auto> typename Super_, auto opts>
class AoSMappedParticles : public AoSParticles<Super_, opts>
{
    using Super  = AoSParticles<Super_, opts>;
    using This   = AoSMappedParticles;
    using Super2 = Super::Super;

    template<template<auto> typename, auto>
    friend class AoSMappedParticles;

public:
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;

    using Particle_t = Super::Particle_t;
    using box_t      = Box<int, dimension>;
    using Super::particles_;
    using Super2::box_;
    using Super2::cellMap_;

    AoSMappedParticles(box_t box = {}, std::size_t size = 0)
        : Super(box, size)
    {
    }

    template<typename Particle_t>
    AoSMappedParticles(box_t box, std::size_t size, Particle_t const& from)
        : Super(box, size, from)
    {
        PHARE_ASSERT(box_.size() > 0);
    }

    template<typename... Args>
    AoSMappedParticles(Args&&... args)
        requires(self_excluding_constructible<This, Super, Args...>())
        : Super{std::forward<Args>(args)...}
    {
    }

    AoSMappedParticles(AoSMappedParticles const&)            = default;
    AoSMappedParticles(AoSMappedParticles&&)                 = default;
    AoSMappedParticles& operator=(AoSMappedParticles&&)      = default;
    AoSMappedParticles& operator=(AoSMappedParticles const&) = default;

    void clear()
    {
        Super::clear();
        cellMap_.clear();
    }

    void erase(box_t const& box)
        requires(storage_mode == StorageMode::VECTOR)
    {
        Super::Super::erase(box);
    }

    template<typename IndexRange>
    auto erase(IndexRange range)
    {
        cellMap_.erase(range);
    }

    template<typename It>
    auto erase(It a, It b)
    {
        return Super::erase(a, b);
    }

    auto& emplace_back()
    {
        auto& part = Super::emplace_back();
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    auto& emplace_back(This const& src, std::size_t const& idx)
    {
        auto& part = Super::emplace_back(src, idx);
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    auto& emplace_back(Particle_t&& p)
    {
        auto& part = Super::emplace_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    auto& emplace_back(Particle_t const& p)
    {
        auto& part = Super::emplace_back(p);
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    template<typename... Args>
    auto& emplace_back(Args const&... args)
    {
        auto& part = Super::emplace_back(args...);
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    template<typename... Args>
    auto& emplace_back(Args&&... args)
    {
        auto& part = Super::emplace_back(args...);
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    void push_back(Particle_t const& p)
    {
        Super::push_back(p);
        cellMap_.add(particles_, particles_.size() - 1);
    }

    void push_back(Particle_t&& p)
    {
        Super::push_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
    }

    void map_particles() const { cellMap_.add(particles_); }
    void empty_map() { cellMap_.empty(); }

    void map_particles(std::size_t const first)
    {
        cellMap_.add(particles_, first, particles_.size());
    }

    void remap()
    {
        empty_map();
        map_particles();
    };

    template<typename... Args>
    auto nbr_particles_in(Args&&... args) const
    {
        return cellMap_.size(args...);
    }

    void export_particles(box_t const& box, auto& dest) const
    {
        cellMap_.export_to(box, *this, dest);
    }

    template<typename Dest, typename Fn>
    void export_particles(box_t const& box, Dest& dest, Fn&& fn) const
    {
        cellMap_.export_to(box, *this, dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        cellMap_.export_to(box, *this, dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(This& dest, Predicate&& pred) const
    {
        cellMap_.export_if(*this, dest, std::forward<Predicate>(pred));
    }

    template<typename Cell>
    void change_icell(Cell const& newCell, std::size_t particleIndex)
    {
        auto oldCell                      = particles_[particleIndex].iCell();
        particles_[particleIndex].iCell() = newCell;
        auto const box_is_valid           = box_.size() > 1;
        if (box_is_valid)
            cellMap_.update(particles_, particleIndex, oldCell);
    }

    template<typename IndexRange, typename Predicate>
    auto partition(IndexRange&& range, Predicate&& pred)
    {
        return cellMap_.partition(range, std::forward<Predicate>(pred));
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cellMap_.print(cell);
    }

    void sortMapping() const { cellMap_.sort(); }

    NO_DISCARD bool is_consistent() const
    {
        if (particles_.size() != cellMap_.size())
            return false;

        for (std::size_t pidx = 0; pidx < particles_.size(); ++pidx)
            if (!cellMap_(particles_[pidx].iCell()).is_indexed(pidx))
                return false;
        return true;
    }

    void check() const;

    // AMR INDEXES!
    NO_DISCARD auto& map(auto const& arr) { return cellMap_(arr); }
    NO_DISCARD auto& map(auto const& arr) const { return cellMap_(arr); }

    auto& box() const { return box_; }
};

template<template<auto> typename Super_, auto opts>
void AoSMappedParticles<Super_, opts>::check() const
{
    if (particles_.size() != cellMap_.size())
    {
        PHARE_LOG_LINE_STR(particles_.size());
        PHARE_LOG_LINE_STR(cellMap_.size());
    }

    core::abort_if(particles_.size() != cellMap_.size());
    core::abort_if(!is_consistent());
}


template<template<auto> typename Super, auto opts>
template<typename T>
struct AoSParticles<Super, opts>::iterator_impl
{
    auto static constexpr dimension = opts.dim;

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using Particle_t        = Super::Particle_t;
    using value_type        = Particle_t;
    using pointer           = Particle_t*;
    using reference         = Particle_t&;

    iterator_impl(T& particles_, std::size_t const s = 0)
        : particles{particles_}
        , curr_pos{s}
    {
    }
    iterator_impl(iterator_impl&& that)      = default;
    iterator_impl(iterator_impl const& that) = default;

    iterator_impl& operator=(iterator_impl&& that)      = default;
    iterator_impl& operator=(iterator_impl const& that) = default;


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
    auto operator+(std::int64_t i) const
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


    auto operator==(iterator_impl const& that) const { return curr_pos == that.curr_pos; }
    auto operator!=(iterator_impl const& that) const { return curr_pos != that.curr_pos; }
    auto operator<(iterator_impl const& that) const { return curr_pos < that.curr_pos; }
    auto operator-(iterator_impl const& that) const { return curr_pos - that.curr_pos; }

    auto& weight() { return deref(particles)[curr_pos].weight(); }
    auto& weight() const { return deref(particles)[curr_pos].weight(); }

    auto& charge() { return deref(particles)[curr_pos].charge(); }
    auto& charge() const { return deref(particles)[curr_pos].charge(); }

    auto& iCell() { return deref(particles)[curr_pos].iCell(); }
    auto& iCell() const { return deref(particles)[curr_pos].iCell(); }

    auto& delta() { return deref(particles)[curr_pos].delta(); }
    auto& delta() const { return deref(particles)[curr_pos].delta(); }

    auto& v() { return deref(particles)[curr_pos].v(); }
    auto& v() const { return deref(particles)[curr_pos].v(); }

    auto& operator*() { return deref(particles)[curr_pos]; }
    auto& operator*() const { return deref(particles)[curr_pos]; }

    auto& operator[](std::size_t const i) { return deref(particles)[curr_pos + i]; }
    auto& operator[](std::size_t const i) const { return deref(particles)[curr_pos + i]; }

    auto copy() const { return deref(particles)[curr_pos]; }

    T particles;
    std::size_t curr_pos;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_AOS_HPP */
