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
        using pointer_tuple_t                  = tuple_fixed_type<ParticleArray*, 3>;

        ParticlesPack(ParticlesPack const& source)            = default;
        ParticlesPack(ParticlesPack&& source)                 = default;
        ParticlesPack& operator=(ParticlesPack&& source)      = default;
        ParticlesPack& operator=(ParticlesPack const& source) = default;

        void setBuffer(ParticlesPack* const source)
        {
            (*this) = source ? *source : ParticlesPack{_name};
        }

        std::string _name;
        ParticleArray* _domainParticles{nullptr};
        ParticleArray* _patchGhostParticles{nullptr};
        ParticleArray* _levelGhostParticles{nullptr};
        ParticleArray* _levelGhostParticlesOld{nullptr};
        ParticleArray* _levelGhostParticlesNew{nullptr};

        auto& name() const { return _name; }

        auto _as_tuple() const
        {
            return std::forward_as_tuple(_domainParticles, _patchGhostParticles,
                                         _levelGhostParticles);
        }

        NO_DISCARD bool isUsable() const
        {
            auto tup = _as_tuple();
            return for_N_all<std::tuple_size_v<pointer_tuple_t>>(
                [&](auto i) { return std::get<i>(tup) != nullptr; });
        }


        NO_DISCARD bool isSettable() const
        {
            auto tup = _as_tuple();
            return for_N_all<std::tuple_size_v<pointer_tuple_t>>(
                [&](auto i) { return std::get<i>(tup) == nullptr; });
        }


        NO_DISCARD auto& domainParticles() const
        {
            if (_domainParticles)
                return *_domainParticles;
            throw std::runtime_error("Error - cannot provide access to domainParticles");
        }
        NO_DISCARD auto& domainParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<ParticlesPack const*>(this)->domainParticles());
        }


        NO_DISCARD auto& patchGhostParticles() const
        {
            if (_patchGhostParticles)
                return *_patchGhostParticles;
            throw std::runtime_error("Error - cannot provide access to patchGhostParticles");
        }
        NO_DISCARD auto& patchGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<ParticlesPack const*>(this)->patchGhostParticles());
        }


        NO_DISCARD auto& levelGhostParticles() const
        {
            if (_levelGhostParticles)
                return *_levelGhostParticles;
            throw std::runtime_error("Error - cannot provide access to levelGhostParticles");
        }
        NO_DISCARD auto& levelGhostParticles()
        {
            return const_cast<ParticleArray&>(
                static_cast<ParticlesPack const*>(this)->levelGhostParticles());
        }


        NO_DISCARD ParticleArray& levelGhostParticlesOld()
        {
            if (_levelGhostParticlesOld)
                return *_levelGhostParticlesOld;
            throw std::runtime_error("Error - cannot provide access to levelGhostParticlesOld");
        }
        NO_DISCARD ParticleArray const& levelGhostParticlesOld() const
        {
            return const_cast<ParticlesPack*>(this)->levelGhostParticlesOld();
        }


        NO_DISCARD ParticleArray& levelGhostParticlesNew()
        {
            if (_levelGhostParticlesNew)
                return *_levelGhostParticlesNew;
            throw std::runtime_error("Error - cannot provide access to levelGhostParticlesNew");
        }
        NO_DISCARD ParticleArray const& levelGhostParticlesNew() const
        {
            return const_cast<ParticlesPack*>(this)->levelGhostParticlesNew();
        }
    };


} // namespace core
} // namespace PHARE

#endif // PARTICLE_PACK_HPP
