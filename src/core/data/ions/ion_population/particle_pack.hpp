#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_PACK_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_PACK_HPP

#include "core/data/particles/particle_array.hpp"

namespace PHARE
{
namespace core
{
    //! ParticlePack is used to pack three ParticleArray together
    /**
     * For convenience, particles are often stored in arrays specific to their
     * location. Ghost particles are particles coming from neighbor patches of
     * the same level, incoming particles are particles coming from coarse-to-fine
     * boundaries, and interior particles are those within the patch physical domain.
     * ParticlesPack conveniently store the three arrays together.
     */
    template<typename ParticleArray>
    struct ParticlesPack
    {
        static constexpr std::size_t dimension = ParticleArray::dimension;
        using particle_array_type              = ParticleArray;

        ParticlesPack(ParticlesPack const& source)            = default;
        ParticlesPack(ParticlesPack&& source)                 = default;
        ParticlesPack& operator=(ParticlesPack&& source)      = default;
        ParticlesPack& operator=(ParticlesPack const& source) = default;

        void setBuffer(ParticlesPack* const source)
        {
            (*this) = source ? *source : ParticlesPack{name_};
        }

        std::string name_;
        ParticleArray* domainParticles{nullptr};
        ParticleArray* patchGhostParticles{nullptr};
        ParticleArray* levelGhostParticles{nullptr};
        ParticleArray* levelGhostParticlesOld{nullptr};
        ParticleArray* levelGhostParticlesNew{nullptr};

        auto& name() const { return name_; }
        NO_DISCARD bool isUsable() const { return domainParticles != nullptr; }
        NO_DISCARD bool isSettable() const { return domainParticles == nullptr; }
    };
} // namespace core
} // namespace PHARE

#endif // PARTICLE_PACK_HPP
