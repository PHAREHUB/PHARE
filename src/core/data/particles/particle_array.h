#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H


#include <cstddef>
#include <utility>
#include <vector>

#include "particle.h"
#include "core/utilities/point/point.h"
#include "core/utilities/cellmap.h"
#include "core/logger.h"

namespace PHARE::core
{
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous              = false;
    static constexpr auto dimension                  = dim;
    static constexpr std::size_t cellmap_bucket_size = 100;
    using Particle_t                                 = Particle<dim>;
    using Vector                                     = std::vector<Particle_t>;
    using iterator                                   = typename Vector::iterator;
    using value_type                                 = Particle_t;
    using box_t                                      = Box<int, dim>;

private:
    using cell_map_t = CellMap<dim, Particle_t, cellmap_bucket_size, int, Point<int, dim>>;

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

    auto begin() const { return particles_.begin(); }
    auto begin()
    {
        clean_ = false;
        return particles_.begin();
    }

    auto end() const { return particles_.end(); }
    auto end()
    {
        clean_ = false;
        return particles_.end();
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
        return particles_.erase(position);
    }
    iterator erase(iterator first, iterator last)
    {
        clean_ = false;
        return particles_.erase(first, last);
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
    Vector particles_;
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

} // namespace core
} // namespace PHARE


#endif
