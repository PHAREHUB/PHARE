#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP


#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/utilities/range/range.hpp"

#include "particle.hpp"


#include <vector>
#include <cstddef>
#include <utility>


namespace PHARE::core
{
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = dim;
    using This                          = ParticleArray<dim>;
    using Particle_t                    = Particle<dim>;
    using Vector                        = std::vector<Particle_t>;

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

    ParticleArray(ParticleArray const& from)            = default;
    ParticleArray(ParticleArray&& from)                 = default;
    ParticleArray& operator=(ParticleArray&& from)      = default;
    ParticleArray& operator=(ParticleArray const& from) = default;

    NO_DISCARD std::size_t size() const { return particles_.size(); }
    NO_DISCARD std::size_t capacity() const { return particles_.capacity(); }

    void clear()
    {
        particles_.clear();
        cellMap_.clear();
    }
    void reserve(std::size_t newSize) { return particles_.reserve(newSize); }
    void resize(std::size_t newSize) { return particles_.resize(newSize); }

    NO_DISCARD auto const& operator[](std::size_t i) const { return particles_[i]; }
    NO_DISCARD auto& operator[](std::size_t i) { return particles_[i]; }

    NO_DISCARD bool operator==(ParticleArray<dim> const& that) const
    {
        return (this->particles_ == that.particles_);
    }

    NO_DISCARD auto begin() const { return particles_.begin(); }
    NO_DISCARD auto begin() { return particles_.begin(); }

    NO_DISCARD auto end() const { return particles_.end(); }
    NO_DISCARD auto end() { return particles_.end(); }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles_.insert(position, first, last);
    }

    NO_DISCARD auto back() { return particles_.back(); }
    NO_DISCARD auto front() { return particles_.front(); }


    auto erase(IndexRange_ range) { cellMap_.erase(range); }

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



    void map_particles() const { cellMap_.add(particles_); }
    void empty_map() { cellMap_.empty(); }


    template<typename Cell>
    bool swap_last_reduce_by_one(Cell const& oldCell, std::size_t const idx)
    {
        // swap last to index
        // swap idx to last
        // --size
        // return true if you need to repeat the current index == expected

        bool const idx_is_last = idx == size() - 1;
        if (!idx_is_last)
        {
            cellMap_.swap(oldCell, particles_[size() - 1].iCell, idx, size() - 1);
            particles_[idx] = particles_[size() - 1];
        }

        cellMap_.erase(*this, oldCell, size() - 1); // doesn't erase from particles vector
        resize(size() - 1);
        return !idx_is_last;
    }

    void swap(std::size_t const a, std::size_t const b)
    {
        cellMap_.swap(particles_[a].iCell, particles_[b].iCell, a, b);
        auto const tmp = particles_[a];
        particles_[a]  = particles_[b];
        particles_[b]  = tmp;
    }


    NO_DISCARD auto nbr_particles_in(box_t const& box) const { return cellMap_.size(box); }

    using cell_t = std::array<int, dim>;
    auto nbr_particles_in(cell_t const& cell) const { return cellMap_.size(cell); }

    void export_particles(box_t const& box, ParticleArray<dim>& dest) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles");
        cellMap_.export_to(box, particles_, dest);
    }

    template<typename Fn>
    void export_particles(box_t const& box, ParticleArray<dim>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (Fn)");
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (box, vector, Fn)");
        cellMap_.export_to(box, particles_.data(), dest, std::forward<Fn>(fn));
    }

    template<typename Predicate>
    void export_particles(ParticleArray& dest, Predicate&& pred) const
    {
        PHARE_LOG_SCOPE(3, "ParticleArray::export_particles (Fn,vector)");
        cellMap_.export_if(particles_.data(), dest, std::forward<Predicate>(pred));
    }


    template<typename Cell>
    void change_icell(Particle_t& /*particle*/, Cell const& oldCell,
                      std::size_t const particleIndex)
    {
        if (!box_.isEmpty())
            cellMap_.update(particles_, particleIndex, oldCell);
    }


    template<typename Cell>
    void change_icell(Cell const& newCell, std::size_t const particleIndex)
    {
        auto const oldCell              = particles_[particleIndex].iCell;
        particles_[particleIndex].iCell = newCell;
        if (!box_.isEmpty())
        {
            cellMap_.update(particles_, particleIndex, oldCell);
        }
    }


    template<typename Predicate>
    auto partition(Predicate&& pred)
    {
        return cellMap_.partition(makeIndexRange(*this), std::forward<Predicate>(pred));
    }

    template<typename Range_t, typename Predicate>
    auto partition(Range_t&& range, Predicate&& pred)
    {
        auto const ret = cellMap_.partition(range, std::forward<Predicate>(pred));
        assert(ret.size() <= range.size());
        return ret;
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cellMap_.print(cell);
    }


    NO_DISCARD bool is_consistent() const
    {
        if (particles_.size() != cellMap_.size())
            return false;

        for (std::size_t pidx = 0; pidx < particles_.size(); ++pidx)
            if (!cellMap_(particles_[pidx].iCell).is_indexed(pidx))
                return false;

        return true;
    }

    void sortMapping() const { cellMap_.sort(); }

    NO_DISCARD auto& vector() { return particles_; }
    NO_DISCARD auto& vector() const { return particles_; }

    auto& box() const { return box_; }


private:
    Vector particles_;
    box_t box_;
    mutable CellMap_t cellMap_;
};

} // namespace PHARE::core


namespace PHARE
{
namespace core
{


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

        NO_DISCARD std::size_t size() const { return weight.size(); }

        template<std::size_t S, typename T>
        NO_DISCARD static std::array<T, S>* _array_cast(T const* array)
        {
            return reinterpret_cast<std::array<T, S>*>(const_cast<T*>(array));
        }

        template<typename Return>
        NO_DISCARD Return _to(std::size_t i)
        {
            return {
                *const_cast<double*>(weight.data() + i),     //
                *const_cast<double*>(charge.data() + i),     //
                *_array_cast<dim>(iCell.data() + (dim * i)), //
                *_array_cast<dim>(delta.data() + (dim * i)), //
                *_array_cast<3>(v.data() + (3 * i)),
            };
        }

        NO_DISCARD auto copy(std::size_t i) { return _to<Particle<dim>>(i); }
        NO_DISCARD auto view(std::size_t i) { return _to<ParticleView<dim>>(i); }

        NO_DISCARD auto operator[](std::size_t i) const { return view(i); }
        NO_DISCARD auto operator[](std::size_t i) { return view(i); }

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

            NO_DISCARD bool operator!=(iterator const& other) const
            {
                return curr_pos != views.size();
            }
            NO_DISCARD auto& operator*() { return views[curr_pos]; }
            NO_DISCARD auto& operator*() const { return views[curr_pos]; }

            std::size_t curr_pos = 0;
            std::vector<ParticleView<dim>> views;
        };

        NO_DISCARD auto as_tuple()
        {
            return std::forward_as_tuple(weight, charge, iCell, delta, v);
        }
        NO_DISCARD auto as_tuple() const
        {
            return std::forward_as_tuple(weight, charge, iCell, delta, v);
        }

        NO_DISCARD auto begin() { return iterator(this); }
        NO_DISCARD auto cbegin() const { return iterator(this); }

        NO_DISCARD auto end() { return iterator(this); }
        NO_DISCARD auto cend() const { return iterator(this); }

        container_t<int> iCell;
        container_t<double> delta;
        container_t<double> weight, charge, v;
    };


    template<std::size_t dim>
    using ContiguousParticlesView = ContiguousParticles<dim, /*OwnedState=*/false>;

} // namespace core
} // namespace PHARE


#endif
