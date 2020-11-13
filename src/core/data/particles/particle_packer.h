#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_H
#define PHARE_CORE_DATA_PARTICLE_PACKER_H


#include <cstddef>
#include <vector>

#include "particle.h"
#include "particle_array.h"

namespace PHARE::core
{
template<std::size_t dim>
class ParticlePacker
{
public:
    ParticlePacker(ParticleArray<dim> const& particles)
        : particles_{particles}
    {
    }

    static auto get(Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    static auto empty()
    {
        Particle<dim> particle;
        return get(particle);
    }

    static auto& keys() { return keys_; }

    auto get(std::size_t i) const { return get(particles_[i]); }
    bool hasNext() const { return it_ < particles_.size(); }
    auto next() { return get(it_++); }

    void pack(ContiguousParticles<dim>& copy)
    {
        auto copyTo = [](auto& a, auto& idx, auto size, auto& v) {
            std::copy(a.begin(), a.begin() + size, v.begin() + (idx * size));
        };
        std::size_t idx = 0;
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
    ParticleArray<dim> const& particles_;
    std::size_t it_ = 0;
    static inline std::array<std::string, 5> keys_{"weight", "charge", "iCell", "delta", "v"};
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
