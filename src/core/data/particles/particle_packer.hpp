#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_HPP
#define PHARE_CORE_DATA_PARTICLE_PACKER_HPP

#include <cstddef>

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class ParticlePacker
{
    using AoSParticles_t = AoSParticles<AoSVector<dim>>;

public:
    ParticlePacker(AoSParticles_t const& particles)
        : particles_{particles}
    {
    }

    static auto get(Particle<dim> const& particle)
    {
        return std::forward_as_tuple(particle.weight_, particle.charge_, particle.iCell_,
                                     particle.delta_, particle.v_);
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

    template<typename SoAParticles_t>
    void pack(SoAParticles_t& copy)
    {
        std::size_t idx = 0;
        while (this->hasNext())
        {
            auto next        = this->next();
            copy.weight(idx) = std::get<0>(next);
            copy.charge(idx) = std::get<1>(next);
            copy.iCell(idx)  = std::get<2>(next);
            copy.delta(idx)  = std::get<3>(next);
            copy.v(idx)      = std::get<4>(next);
            idx++;
        }
    }

private:
    AoSParticles_t const& particles_;
    std::size_t it_ = 0;
    static inline std::array<std::string, 5> keys_{"weight", "charge", "iCell", "delta", "v"};
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
