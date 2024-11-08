#ifndef PHARE_PARTICLE_TILE_HPP
#define PHARE_PARTICLE_TILE_HPP

#include "core/utilities/box/box.hpp"
#include "core/data/clusters/clusters.hpp"

namespace PHARE::core
{

template<typename ParticleArray>
struct ParticleCluster : public Box<int, ParticleArray::dimension>
{
    // *this: box constructed first
    ParticleArray domainParticles_{*this};
    ParticleArray patchGhostParticles_{*this};
    ParticleArray levelGhostParticles_{*this};
    ParticleArray levelGhostParticlesNew_{*this};
    ParticleArray levelGhostParticlesOld_{*this};
};


template<typename ParticleArray>
class ParticleClusterSet : public ClusterSet<ParticleCluster<ParticleArray>>
{
public:
    static auto constexpr dimension = ParticleArray::dimension;

    auto select_particles(Box<int, dimension> const& from) const
    {
        ParticleArray selected;

        auto overlaped_clusters = export_overlaped_with(from);
        for (auto const& [complete, cluster] : overlaped_clusters)
        {
            if (complete)
            {
                selected.push_back(cluster.domainParticles_);
            }
            else
            {
                auto intersection = cluster * from;
                cluster.domainParticles.export_to(*intersection, selected);
            }
        }
        return selected;
    }
};
} // namespace PHARE::core

#endif
