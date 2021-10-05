#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H

#include <vector>
#include <cstddef>
#include "core/data/particles/particle.h"

namespace PHARE::core
{
template<std::size_t dim>
class ParticleArray
{
public:
    static constexpr bool is_contiguous = false;
    static constexpr auto dimension     = dim;
    using Particle_t                    = core::Particle<dim>;
    using This                          = ParticleArray<dim>;
    using Vector                        = std::vector<Particle_t>;
    using iterator                      = typename Vector::iterator;
    using value_type                    = Particle_t;

    ParticleArray() {}
    ParticleArray(std::size_t size)
        : particles(size)
    {
    }

    ParticleArray(ParticleArray const& that)
        : particles(that.particles)
    {
    }

    ParticleArray(ParticleArray&& that)
        : particles(std::move(that.particles))
    {
    }

    ParticleArray(std::size_t size, Particle_t&& particle)
        : particles(size, particle)
    {
    }

    auto& operator=(ParticleArray const& that)
    {
        this->particles = that.particles;
        return *this;
    }
    auto& operator=(ParticleArray&& that)
    {
        this->particles = std::move(that.particles);
        return *this;
    }

    auto& operator=(Vector&& vector)
    {
        this->particles = std::move(vector);
        return *this;
    }

    std::size_t size() const { return particles.size(); }

    void clear() { return particles.clear(); }
    void reserve(std::size_t newSize) { return particles.reserve(newSize); }
    void resize(std::size_t newSize) { return particles.resize(newSize); }

    auto& operator[](std::size_t i) const { return particles[i]; }
    auto& operator[](std::size_t i) { return particles[i]; }

    bool operator==(This const& that) const { return (this->particles == that.particles); }

    auto begin() const { return particles.begin(); }
    auto begin() { return particles.begin(); }

    auto end() const { return particles.end(); }
    auto end() { return particles.end(); }

    template<class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    {
        particles.insert(position, first, last);
    }

    auto data() { return particles.data(); }
    auto data() const { return particles.data(); }

    auto back() { return particles.back(); }
    auto back() const { return particles.back(); }

    auto front() { return particles.front(); }
    auto front() const { return particles.front(); }

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

template<std::size_t dim>
void append(ParticleArray<dim> const& src, ParticleArray<dim>& dst)
{
    std::copy(std::begin(src), std::end(src), std::back_inserter(dst));
}

} // namespace PHARE::core


#endif
