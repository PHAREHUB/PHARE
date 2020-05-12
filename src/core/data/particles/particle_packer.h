#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_H
#define PHARE_CORE_DATA_PARTICLE_PACKER_H


#include <cstddef>
#include <vector>

#include "particle.h"
#include "particle_array.h"

namespace PHARE::core
{
template<size_t dim, typename Float = double>
class ParticlePacker
{
    using Particle_t            = Particle<dim, Float>;
    using ParticleArray_t       = ParticleArray<dim, Float>;
    using ContiguousParticles_t = ContiguousParticles<dim, Float>;

public:
    ParticlePacker(ParticleArray_t const& particles)
        : particles_{particles}
    {
    }

    static auto get(Particle_t const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    static auto empty()
    {
        Particle_t particle;
        return get(particle);
    }

    static auto& keys() { return keys_; }

    auto get(size_t i) const { return get(particles_[i]); }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

    void pack(ContiguousParticles_t& copy)
    {
        auto copyTo = [](auto& a, auto& idx, auto size, auto& v) {
            std::copy(a.begin(), a.begin() + size, v.begin() + (idx * size));
        };
        size_t idx = 0;
        while (this->hasNext())
        {
            auto next        = this->next();
            copy.weight[idx] = std::get<0>(next);
            copy.charge[idx] = std::get<1>(next);
            copyTo(std::get<2>(next), idx, dim, copy.iCell);
            copyTo(std::get<3>(next), idx, dim, copy.delta);
            copyTo(std::get<4>(next), idx, 3, copy.v);
            idx++;
        }
    }

private:
    ParticleArray_t const& particles_;
    size_t it_ = 0;
    static inline std::array<std::string, 5> keys_{"weight", "charge", "iCell", "delta", "v"};
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
