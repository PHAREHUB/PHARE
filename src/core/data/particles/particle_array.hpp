#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H


#include <cstddef>
#include <utility>
#include <vector>

#include "core/utilities/bucketlist.hpp"
#include "particle.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/cellmap.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/iterators.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous              = false;
    static constexpr auto dimension                  = dim;
    static constexpr std::size_t cellmap_bucket_size = 100;
    using This                                       = ParticleArray<dim>;
    using Particle_t                                 = Particle<dim, BucketListItem>;
    using Vector                                     = std::vector<Particle_t>;

private:
    using cell_map_t = CellMap<dim, Particle_t, cellmap_bucket_size, int, Point<int, dim>>;

    template<typename T>
    struct ParticleArrayIterator : wrapped_iterator<T, Vector>
    {
        using Super      = wrapped_iterator<T, Vector>;
        using iterator_t = typename Super::iterator_t;

        ParticleArrayIterator(iterator_t iter, T& container, cell_map_t& cellmap)
            : wrapped_iterator<T, Vector>{iter, &container}
            , cm{&cellmap}
        {
        }

        template<typename Cell>
        void change_icell(Cell const& newCell)
        {
            // cm->update(*(this->it), newCell);
            this->it->iCell = newCell;
        }

        bool operator==(ParticleArrayIterator const& other) const
        {
            bool superbool = (static_cast<Super const&>(*this) == static_cast<Super const&>(other));
            bool cmb       = (cm == other.cm);
            return superbool and cmb;
        }
        ParticleArrayIterator& operator=(ParticleArrayIterator const& other) = default;
        cell_map_t* cm;
    };

public:
    using value_type     = Particle_t;
    using box_t          = Box<int, dim>;
    using iterator       = ParticleArrayIterator<This>;
    using const_iterator = ParticleArrayIterator<This const>;


public:
    ParticleArray() {}
    ParticleArray(std::size_t size)
        : particles_(size)
    {
    }

    ParticleArray(std::size_t size, Particle_t&& particle)
        : particles_(size, particle)
    {
    }

    std::size_t size() const { return particles_.size(); }
    std::size_t capacity() const { return particles_.capacity(); }

    void clear()
    {
        clean_ = false;
        return particles_.clear();
    }
    void reserve(std::size_t newSize)
    {
        clean_ = false;
        return particles_.reserve(newSize);
    }
    void resize(std::size_t newSize)
    {
        clean_ = false;
        return particles_.resize(newSize);
    }

    auto const& operator[](std::size_t i) const { return particles_[i]; }

    bool operator==(ParticleArray<dim> const& that) const
    {
        return (this->particles_ == that.particles_);
    }

    auto begin() const { return const_iterator{particles_.begin(), *this, cell_map_}; }
    auto begin()
    {
        clean_ = false;
        return iterator{particles_.begin(), *this, cell_map_};
    }

    auto end() const { return const_iterator{particles_.end(), *this, cell_map_}; }
    auto end()
    {
        clean_ = false;
        return iterator{particles_.end(), *this, cell_map_};
    }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        clean_ = false;
        particles_.insert(position, first, last);
    }

    auto back()
    {
        clean_ = false;
        return particles_.back();
    }
    auto front()
    {
        clean_ = false;
        return particles_.front();
    }

    iterator erase(iterator position)
    {
        clean_ = false;
        return iterator{particles_.erase(position), *this, cell_map_};
    }
    iterator erase(iterator first, iterator last)
    {
        clean_ = false;
        return iterator{particles_.erase(first.it, last.it), *this, cell_map_};
    }

    Particle_t& emplace_back()
    {
        clean_ = false;
        return particles_.emplace_back();
    }
    Particle_t& emplace_back(Particle_t&& p)
    {
        clean_ = false;
        return particles_.emplace_back(p);
    }
    Particle_t& emplace_back(Particle_t const& p)
    {
        clean_ = false;
        return particles_.emplace_back(p);
    }

    void push_back(Particle_t const& p)
    {
        clean_ = false;
        particles_.push_back(p);
    }
    void push_back(Particle_t&& p)
    {
        clean_ = false;
        particles_.push_back(std::forward<Particle_t>(p));
    }

    void swap(ParticleArray<dim>& that)
    {
        clean_ = false;
        std::swap(this->particles_, that.particles_);
    }

    void map_particles() const
    {
        cell_map_.add(particles_);
        clean_ = true;
    }
    void empty_map() { cell_map_.empty(); }


    auto nbr_particles_in(box_t const& box) const
    {
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        auto s = cell_map_.size(box);
        return s;
    }


    auto select(box_t const& box) const
    {
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        return cell_map_.select(box);
    }


    void export_particles(box_t const& box, ParticleArray<dim>& dest) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles");
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        cell_map_.export_to(box, dest.particles_);
    }

    template<typename Fn>
    void export_particles(box_t const& box, ParticleArray<dim>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn)");
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        cell_map_.export_to(box, dest.particles_, std::forward<Fn>(fn));
    }

    template<typename Fn>
    void export_particles(box_t const& box, std::vector<Particle_t>& dest, Fn&& fn) const
    {
        PHARE_LOG_SCOPE("ParticleArray::export_particles (Fn,vector)");
        if (!clean_)
        {
            cell_map_.empty();
            map_particles();
        }
        cell_map_.export_to(box, dest, std::forward<Fn>(fn));
    }

    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        cell_map_.print(cell);
    }

private:
    bool mutable clean_{false};
    mutable Vector particles_;
    mutable cell_map_t cell_map_;
};

} // namespace PHARE::core


namespace PHARE
{
namespace core
{
    template<std::size_t dim>
    void empty(ParticleArray<dim>& array)
    {
        array.clear();
    }

    template<std::size_t dim>
    void swap(ParticleArray<dim>& array1, ParticleArray<dim>& array2)
    {
        array1.swap(array2);
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

        auto copy(std::size_t i) { return _to<Particle<dim, BucketListItem>>(i); }
        auto view(std::size_t i) { return _to<ParticleView<dim, BucketListItem>>(i); }

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
            std::vector<ParticleView<dim, BucketListItem>> views;
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

} // namespace core
} // namespace PHARE


#endif
