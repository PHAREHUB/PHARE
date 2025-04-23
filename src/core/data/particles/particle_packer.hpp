#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_HPP
#define PHARE_CORE_DATA_PARTICLE_PACKER_HPP


#include "core/def.hpp"
#include "particle.hpp"
#include "particle_array.hpp"

#include <cstddef>

namespace PHARE::core
{
template<std::size_t dim>
class ParticlePacker
{
    constexpr static Particle<dim> default_particle{};

public:
    static constexpr std::size_t n_keys = 5;

    ParticlePacker(ParticleArray<dim> const& particles)
        : particles_{particles}
    {
    }

    NO_DISCARD static constexpr auto get(Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.charge, particle.iCell,
                                     particle.delta, particle.v);
    }

    static constexpr auto empty() { return get(default_particle); }

    // sometimes we use this to infer the size of an ParticleArray
    // could be "charge" either
    NO_DISCARD static auto arbitrarySingleValueKey() { return "weight"; }

    NO_DISCARD static auto& keys() { return keys_; }

    NO_DISCARD auto get(std::size_t i) const { return get(particles_[i]); }
    NO_DISCARD bool hasNext() const { return it_ < particles_.size(); }
    NO_DISCARD auto next() { return get(it_++); }

    void pack(ContiguousParticles<dim>& soa) const
    {
        for (auto const& particle : particles_)
            soa.push_back(particle);
    }

    template<typename Fn>
    void pack_ranges_into(Fn const fn, std::size_t const S = 2048) const
    {
        ContiguousParticles<dim> soa{S};
        soa.clear(); // reserved but 0 size

        std::size_t i = 0;
        for (; i < particles_.size() / S; ++i)
        {
            std::size_t const pi = i * S;
            for (std::size_t bi = 0; bi < S; ++bi)
                soa.push_back(particles_[pi + bi]);
            fn(soa, pi);
            soa.clear();
        }

        if (auto const remaining = particles_.size() % S)
        {
            std::size_t const pi = i * S;
            for (std::size_t bi = 0; bi < remaining; ++bi)
                soa.push_back(particles_[pi + bi]);
            fn(soa, pi);
        }
    }

private:
    ParticleArray<dim> const& particles_;
    std::size_t it_ = 0;
    static inline std::array<std::string, n_keys> const keys_{"weight", "charge", "iCell", "delta",
                                                              "v"};
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
