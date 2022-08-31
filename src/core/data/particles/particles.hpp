#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include <vector>
#include <cstddef>
#include <utility>

#include "particle.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/utilities/indexer.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/range/range.hpp"

#include "core/def/types.hpp"


namespace PHARE::core
{
template<typename Particle, typename Allocator = typename std::vector<Particle>::allocator_type>
class ParticleArray
{
public:
    static constexpr bool is_host_mem   = true;
    static constexpr bool is_contiguous = false;
    static constexpr auto dim           = Particle::dimension;
    static constexpr auto dimension     = Particle::dimension;
    using This                          = ParticleArray<Particle, Allocator>;
    using Particle_t                    = Particle;
    using Vector                        = std::vector<Particle_t, Allocator>;

private:
    using CellMap_t   = CellMap<dim, int>;
    using IndexRange_ = IndexRange<This>;


public:
    using value_type     = Particle_t;
    using box_t          = Box<int, dim>;
    using iterator       = typename Vector::iterator;
    using const_iterator = typename Vector::const_iterator;



public:
    ParticleArray(box_t box)
        : box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
    }

    ParticleArray(box_t box, std::size_t size)
        : particles_(size)
        , box_{box}
        , cellMap_{box_}
    {
        assert(box_.size() > 0);
    }

    ParticleArray(ParticleArray const& from) = default;
    ParticleArray(ParticleArray&& from)      = default;
    ParticleArray& operator=(ParticleArray&& from) = default;
    ParticleArray& operator=(ParticleArray const& from) = default;

    std::size_t size() const { return particles_.size(); }
    std::size_t capacity() const { return particles_.capacity(); }

    // template<typename Allocator0>
    // ParticleArray(ParticleArray<Particle, Allocator0> const& that)
    //     : ParticleArray{0}
    // {
    //     PHARE::Vector<Particle_t>::copy(this->particles_, that.particles_);
    // }


    // template<typename Allocator0>
    // ParticleArray(ParticleArray<Particle, Allocator0>&& that)
    //     : ParticleArray{0}
    // {
    //     PHARE::Vector<Particle_t>::move(this->particles_, that.particles_);
    // }

    // template<typename Allocator0>
    // auto& operator=(ParticleArray<Particle, Allocator0> const& that)
    // {
    //     PHARE::Vector<Particle_t>::copy(this->particles_, that.particles_);
    //     return *this;
    // }


    // template<typename Allocator0>
    // auto& operator=(ParticleArray<Particle, Allocator0>&& that)
    // {
    //     PHARE::Vector<Particle_t>::move(this->particles_, that.particles_);
    //     return *this;
    // }

    // template<typename Allocator0>
    // auto& operator=(std::vector<Particle, Allocator0>&& particles)
    // {
    //     PHARE::Vector<Particle_t>::move(this->particles_, particles);
    //     cellMap_.clear();
    //     map_particles();
    //     return *this;
    // }

    void clear()
    {
        particles_.clear();
        cellMap_.clear();
    }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }

    auto const& operator[](std::size_t i) const { return particles_[i]; }
    auto& operator[](std::size_t i) { return particles_[i]; }


    bool operator==(This const& that) const { return (this->particles_ == that.particles_); }

    auto begin() const { return particles_.begin(); }
    auto begin() { return particles_.begin(); }

    auto end() const { return particles_.end(); }
    auto end() { return particles_.end(); }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles_.insert(position, first, last);
    }

    auto& data() { return particles_.data(); }
    auto& data() const { return particles_.data(); }

    auto& back() { return particles_.back(); }
    auto& back() const { return particles_.back(); }

    auto& front() { return particles_.front(); }
    auto& front() const { return particles_.front(); }

    auto erase(IndexRange_& range) { cellMap_.erase(particles_, range); }
    auto erase(IndexRange_&& range)
    {
        // TODO move ctor for range?
        cellMap_.erase(std::forward<IndexRange_>(range));
    }

    iterator erase(iterator first, iterator last)
    {
        // should we erase particles indexes associated with these iterators from the cellmap?
        // probably it does not matter if not. The reason is that
        // particles erased from the particlearray are so because they left
        // the patch cells to an outside cell.
        // But in principle that cell will never be accessed because it is outside the patch.
        // The only thing "bad" if these indexes are not deleted is that the
        // size of the cellmap becomes unequal to the size of the particleArray.
        // but  ¯\_(ツ)_/¯
        return particles_.erase(first, last);
    }


    Particle_t& emplace_back()
    {
        auto& part = particles_.emplace_back();
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }



    Particle_t& emplace_back(Particle_t&& p)
    {
        auto& part = particles_.emplace_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
        return part;
    }

    void push_back(Particle_t const& p)
    {
        particles_.push_back(p);
        cellMap_.add(particles_, particles_.size() - 1);
    }

    void push_back(Particle_t&& p)
    {
        particles_.push_back(std::forward<Particle_t>(p));
        cellMap_.add(particles_, particles_.size() - 1);
    }

    // template<typename... Args>
    // Particle_t& emplace_back(Args const&... args)
    // {
    //     auto& part = particles_.emplace_back(args...);
    //     cellMap_.add(particles_, particles_.size() - 1);
    //     return part;
    // }
    // template<typename... Args>
    // Particle_t& emplace_back(Args&&... args)
    // {
    //     auto& part = particles_.emplace_back(Particle_t{args...});
    //     cellMap_.add(particles_, particles_.size() - 1);
    //     return part;
    // }

    void swap(This& that) { std::swap(this->particles_, that.particles_); }

    void map_particles() const { cellMap_.add(particles_); }
    void empty_map() { cellMap_.empty(); }


    auto nbr_particles_in(box_t const& box) const { return cellMap_.size(box); }

    void export_particles(box_t const& box, This& dest) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles");
        cellMap_.export_to(box, particles_, dest);
    }

    template<typename Fn>
    void export_particles(box_t const& box, This& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn)");
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (box, vector, Fn)");
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(This& dest, Predicate&& pred) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn,vector)");
        cellMap_.export_if(particles_.data(), dest, std::forward<Predicate>(pred));
    }


    template<typename Cell>
    void change_icell(Cell const& newCell, std::size_t particleIndex)
    {
        auto oldCell                    = particles_[particleIndex].iCell;
        particles_[particleIndex].iCell = newCell;
        if (!box_.isEmpty())
        {
            cellMap_.update(particles_, particleIndex, oldCell);
        }
    }


    template<typename Predicate>
    auto partition(Predicate&& pred)
    {
        PHARE_LOG_LINE_STR(particles_.size());
        assert(particles_.size() > 0);

        check();
        return cellMap_.partition(makeIndexRange(*this), std::forward<Predicate>(pred));
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cellMap_.print(cell);
    }


    bool is_mapped() const
    {
        bool ok = true;

        check();

        for (std::size_t pidx = 0; pidx < particles_.size(); ++pidx)
        {
            auto const& p = particles_[pidx];
            auto& icell   = p.iCell;
            auto l        = cellMap_.list_at(icell);
            if (!l)
                throw std::runtime_error("particle cell not mapped");
            auto& ll = l->get();
            if (!ll.is_indexed(pidx))
                throw std::runtime_error("particle not indexed");
        }
        return true;
    }

    void sortMapping() const { cellMap_.sort(); }

    auto& vector() { return particles_; }
    auto& vector() const { return particles_; }

    auto& box() const { return box_; }


    void check()
    {
        PHARE_LOG_LINE_STR(particles_.size());
        PHARE_LOG_LINE_STR(cellMap_.size());

        core::abort_if(particles_.size() != cellMap_.size());
        // throw std::runtime_error("particle array not mapped, map.size() != array.size()");
    }

private:
    Vector particles_;
    box_t box_;
    mutable CellMap_t cellMap_;

    // allow other versions of this class to have access to each others particle vector
    template<typename, typename>
    friend class ParticleArray;
};


template<typename Particle, typename Allocator>
void empty(ParticleArray<Particle, Allocator>& array)
{
    array.clear();
}

template<typename Particle, typename Alloc0, typename Alloc1>
void swap(ParticleArray<Particle, Alloc0>& array1, ParticleArray<Particle, Alloc1>& array2)
{
    static_assert(std::is_same_v<Alloc0, Alloc1>);

    array1.swap(array2);
}

template<typename Particle, typename Allocator>
void append(ParticleArray<Particle, Allocator> const& src, ParticleArray<Particle, Allocator>& dst)
{
    std::copy(std::begin(src), std::end(src), std::back_inserter(dst));
}


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

    auto as_tuple() { return std::forward_as_tuple(weight, charge, iCell, delta, v); }
    auto as_tuple() const { return std::forward_as_tuple(weight, charge, iCell, delta, v); }

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

} // namespace PHARE::core


#endif
