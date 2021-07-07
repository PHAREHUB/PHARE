#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H


#include <cstddef>
#include <vector>

#include "particle.h"

namespace PHARE::core
{
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = dim;
    using Particle_t                    = Particle<dim>;
    using Vector                        = std::vector<Particle_t>;
    using iterator                      = typename Vector::iterator;
    using value_type                    = Particle_t;

    ParticleArray() {}
    ParticleArray(std::size_t size)
        : particles(size)
    {
    }

    ParticleArray(std::size_t size, Particle_t&& particle)
        : particles(size, particle)
    {
    }

    std::size_t size() const { return particles.size(); }

    void clear() { return particles.clear(); }
    void reserve(std::size_t newSize) { return particles.reserve(newSize); }
    void resize(std::size_t newSize) { return particles.resize(newSize); }

    auto& operator[](std::size_t i) const { return particles[i]; }
    auto& operator[](std::size_t i) { return particles[i]; }

    bool operator==(ParticleArray<dim> const& that) const
    {
        return (this->particles == that.particles);
    }

    auto begin() const { return particles.begin(); }
    auto begin() { return particles.begin(); }

    auto end() const { return particles.end(); }
    auto end() { return particles.end(); }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles.insert(position, first, last);
    }

    auto back() { return particles.back(); }
    auto front() { return particles.front(); }

    iterator erase(iterator position) { return particles.erase(position); }
    iterator erase(iterator first, iterator last) { return particles.erase(first, last); }

    Particle_t& emplace_back() { return particles.emplace_back(); }
    Particle_t& emplace_back(Particle_t&& p) { return particles.emplace_back(p); }

    void push_back(Particle_t const& p) { particles.push_back(p); }
    void push_back(Particle_t&& p) { particles.push_back(p); }

    void swap(ParticleArray<dim>& that) { std::swap(this->particles, that.particles); }

private:
    Vector particles;
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
