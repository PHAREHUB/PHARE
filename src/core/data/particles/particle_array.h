#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H


#include <cstddef>
#include <vector>

#include "particle.h"

namespace PHARE::core
{
template<std::size_t dim, bool contiguous = false>
class ParticleArray;

template<std::size_t dim>
struct ParticleArray<dim, false>
{
    using value_type                   = Particle<dim, false>;
    using iterator                     = typename std::vector<value_type>::iterator;
    static constexpr bool is_contigous = false;

    ParticleArray() {}

    ParticleArray(size_t size)
        : particles{size}
    {
    }

    auto& operator[](size_t idx) { return particles[idx]; }
    auto& operator[](size_t idx) const { return particles[idx]; }

    template<typename T>
    auto& emplace_back()
    {
        return particles.emplace_back();
    }

    template<typename... Ts>
    auto emplace_back(Ts&&... args)
    {
        return particles.emplace_back(std::forward<Ts>(args)...);
    }

    template<typename T>
    void push_back(T t)
    {
        particles.push_back(t);
    }

    auto begin() const { return particles.begin(); }
    auto begin() { return particles.begin(); }
    auto end() const { return particles.end(); }
    auto end() { return particles.end(); }

    void clear() { particles.clear(); }

    template<typename A, typename B>
    void erase(A a, B b)
    {
        particles.erase(a, b);
    }

    template<typename A, typename B, typename C>
    void insert(A a, B b, C c)
    {
        particles.insert(a, b, c);
    }

    auto size() const { return particles.size(); }

    std::vector<value_type> particles;

    ParticleArray(const ParticleArray&)             = delete;
    ParticleArray(const ParticleArray&&)            = delete;
    ParticleArray& operator&(const ParticleArray&)  = delete;
    ParticleArray& operator&(const ParticleArray&&) = delete;
};

template<std::size_t dim> // unfinished
struct ParticleArray<dim, true>
{
    static constexpr bool is_contigous = true;

    template<typename T>
    void emplace_back(T t)
    {
    }
    template<typename T>
    void push_back(T t)
    {
    }

    auto begin() const { return 0; }
    auto begin() { return 0; }
    auto end() const { return 0; }
    auto end() { return 0; }

    void clear()
    {
        std::apply(particles.vectorTuple(), [](auto& arg) { arg.clear(); });
        particles.size_ = particles.idx_ = 0;
    }

    Particle<dim, true> particles;

    ParticleArray(const ParticleArray&)             = delete;
    ParticleArray(const ParticleArray&&)            = delete;
    ParticleArray& operator&(const ParticleArray&)  = delete;
    ParticleArray& operator&(const ParticleArray&&) = delete;
};

template<std::size_t dim, bool contiguous = false>
void empty(ParticleArray<dim, contiguous>& array)
{
    array.clear();
}

template<std::size_t dim, bool contiguous = false>
void swap(ParticleArray<dim, contiguous>& array1, ParticleArray<dim, contiguous>& array2)
{
    if constexpr (!contiguous)
        std::swap(array1.particles, array2.particles);
}

} // namespace PHARE::core


#endif
