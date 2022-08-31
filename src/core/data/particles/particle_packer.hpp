#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_HPP
#define PHARE_CORE_DATA_PARTICLE_PACKER_HPP

#include <vector>
#include <cstddef>

#include "particle.hpp"
#include "particle_array.hpp"
#if defined(PHARE_HAVE_UMPIRE)
// #include "core/data/particles/llnl/particle_array.h"
#include "core/def/types.hpp"
#endif

namespace PHARE::core
{
// PGI compiler (nvc++ 21.3-0) doesn't like static initializations of arrays,
//   would result in empty strings
inline std::array<std::string, 5> packer_keys()
{
    // The order of this array must match the tuple order of ParticlePacker::get(particle)
    return {"weight", "charge", "iCell", "delta", "v"};
}

template<std::size_t dim>
class ParticlePacker
{
    using Particle_t = Particle<dim>;

public:
    template<typename Allocator>
    ParticlePacker(ParticleArray<Particle_t, Allocator> const& particles)
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


    // sometimes we use this to infer the size of an ParticleArray
    // could be "charge" either
    static auto arbitrarySingleValueKey() { return "weight"; }

    static auto keys() { return packer_keys(); }

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
#if defined(PHARE_HAVE_UMPIRE)
    ParticleArray<Particle_t> particles_;
#else
    ParticleArray<Particle_t> const& particles_;
#endif
    std::size_t it_ = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
