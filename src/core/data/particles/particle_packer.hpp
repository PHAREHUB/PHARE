#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_HPP
#define PHARE_CORE_DATA_PARTICLE_PACKER_HPP


#include <cstddef>
#include <vector>

#include "particle.hpp"
#include "particle_array.hpp"
#include "core/def.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class ParticlePacker
{
    constexpr static Particle<dim> default_particle{};

public:
    static constexpr std::size_t n_keys = 4;

    ParticlePacker(ParticleArray<dim> const& particles)
        : particles_{particles}
    {
    }

    NO_DISCARD static auto get(Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight, particle.iCell, particle.delta, particle.v);
    }

    static constexpr auto empty() { return get(default_particle); }

    // sometimes we use this to infer the size of an ParticleArray
    NO_DISCARD static auto arbitrarySingleValueKey() { return "weight"; }

    NO_DISCARD static auto& keys() { return keys_; }

    NO_DISCARD auto get(std::size_t i) const { return get(particles_[i]); }
    NO_DISCARD bool hasNext() const { return it_ < particles_.size(); }
    NO_DISCARD auto next() { return get(it_++); }

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
            copyTo(std::get<1>(next), idx, dim, copy.iCell);
            copyTo(std::get<2>(next), idx, dim, copy.delta);
            copyTo(std::get<3>(next), idx, 3, copy.v);
            idx++;
        }
    }

private:
    ParticleArray<dim> const& particles_;
    std::size_t it_ = 0;
    static inline std::array<std::string, n_keys> const keys_{"weight", "iCell", "delta", "v"};

    static_assert(std::tuple_size_v<decltype(get(default_particle))> == n_keys);
    static_assert(std::tuple_size_v<decltype(std::declval<ContiguousParticles<dim>&>().as_tuple())>
                  == n_keys);
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
