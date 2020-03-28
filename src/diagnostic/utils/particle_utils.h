#ifndef PARTICLE_UTILS_H
#define PARTICLE_UTILS_H


#include "core/data/particles/particle_array.h"

namespace PHARE::diagnostic {



template<size_t dim>
class ParticlePacker
{
public:
    ParticlePacker(core::ParticleArray<dim> const& particles)
        : particles_{particles}
    {
    }

    static auto get(core::Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    static auto empty()
    {
        core::Particle<dim> particle;
        return get(particle);
    }

    static auto& keys() { return keys_; }

    auto get(size_t i) const { return get(particles_[i]); }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

private:
    core::ParticleArray<dim> const& particles_;
    size_t it_ = 0;
    static inline std::array<std::string, 5> keys_{"weight", "charge", "iCell", "delta", "v"};
};



template<std::size_t dim>
struct ContiguousParticles
{
    std::vector<int> iCell;
    std::vector<float> delta;
    std::vector<double> weight, charge, v;
    size_t size_;
    auto size() const { return size_; }
    ContiguousParticles(size_t s)
        : iCell(s * dim)
        , delta(s * dim)
        , weight(s)
        , charge(s)
        , v(s * 3)
        , size_(s)
    {
    }
};

}
#endif // PARTICLE_UTILS_H
