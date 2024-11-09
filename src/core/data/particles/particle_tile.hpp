#ifndef PHARE_PARTICLE_TILE_HPP
#define PHARE_PARTICLE_TILE_HPP

#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"

namespace PHARE::core
{

template<typename ParticleArray>
struct ParticleTile : public Box<int, ParticleArray::dimension>
{
    // *this: box constructed first
    ParticleArray domainParticles_{*this};
    ParticleArray patchGhostParticles_{*this};
    ParticleArray levelGhostParticles_{*this};
    ParticleArray levelGhostParticlesNew_{*this};
    ParticleArray levelGhostParticlesOld_{*this};
};


template<typename ParticleArray>
class ParticleTileSet : public TileSet<ParticleTile<ParticleArray>>
{
public:
    static auto constexpr dimension = ParticleArray::dimension;

    auto select_particles(Box<int, dimension> const& from) const
    {
        ParticleArray selected;

        auto overlaped_tiles = export_overlaped_with(from);
        for (auto const& [complete, tile] : overlaped_tiles)
        {
            if (complete)
            {
                selected.push_back(tile.domainParticles_);
            }
            else
            {
                auto intersection = tile * from;
                tile.domainParticles.export_to(*intersection, selected);
            }
        }
        return selected;
    }
};
} // namespace PHARE::core

#endif
