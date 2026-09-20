#ifndef PHARE_CORE_DATA_PARTICLE_PACKER_HPP
#define PHARE_CORE_DATA_PARTICLE_PACKER_HPP

#include <cstddef>

#include "core/def.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{
// PGI compiler (nvc++ 21.3-0) doesn't like static initializations of arrays,
//   would result in empty strings
inline std::array<std::string, 5> packer_keys()
{
    // The order of this array must match the tuple order of ParticlePacker::get(particle)
    return {"weight", "charge", "iCell", "delta", "v"};
}

template<typename SoAParticles_t, typename Particle_t>
void pack_one(SoAParticles_t& copy, std::size_t const idx, Particle_t const& particle)
{
    copy.weight(idx) = particle.weight();
    copy.charge(idx) = particle.charge();
    copy.iCell(idx)  = particle.iCell();
    copy.delta(idx)  = particle.delta();
    copy.v(idx)      = particle.v();
}

// per-layout size()/pack() bodies, specialized on LayoutMode instead of branching
template<LayoutMode layout_mode, typename ParticleArray_>
struct PackerBackend
{
    template<ParticleType>
    NO_DISCARD static std::size_t size(ParticleArray_ const& particles)
    {
        return particles.size();
    }

    template<ParticleType, typename SoAParticles_t>
    static void pack(ParticleArray_ const& particles, SoAParticles_t& copy)
    {
        for (std::size_t idx = 0; idx < particles.size(); ++idx)
            pack_one(copy, idx, particles[idx]);
    }
};

template<typename ParticleArray_>
struct PackerBackend<LayoutMode::SoA, ParticleArray_>
{
    template<ParticleType>
    NO_DISCARD static std::size_t size(ParticleArray_ const& particles)
    {
        return particles.size();
    }

    template<ParticleType, typename SoAParticles_t>
    static void pack(ParticleArray_ const& particles, SoAParticles_t& copy)
    {
        copy = particles;
    }
};

template<typename ParticleArray_>
struct PackerBackend<LayoutMode::AoSPCTS, ParticleArray_>
{
    // adjacent tiles duplicate cells in their private ghost halos, so LevelGhost
    // counts/packs only the clamp-owner's copy (TileSet::tag_cells_) to match
    // AoSMapped. Domain has no such duplication
    template<ParticleType ptype>
    NO_DISCARD static std::size_t size(ParticleArray_ const& particles)
    {
        if constexpr (ptype != ParticleType::LevelGhost)
            return particles.size();
        else
        {
            std::size_t n = 0;
            visit<ptype>(particles, [&](auto const& pc, auto const& amr) {
                n += pc(pc.local_cell(amr)).size();
            });
            return n;
        }
    }

    template<ParticleType ptype, typename SoAParticles_t>
    static void pack(ParticleArray_ const& particles, SoAParticles_t& copy)
    {
        std::size_t idx = 0;
        visit<ptype>(particles, [&](auto const& pc, auto const& amr) {
            for (auto const& particle : pc(pc.local_cell(amr)))
                pack_one(copy, idx++, particle);
        });
        if constexpr (ptype != ParticleType::LevelGhost)
            assert(idx == particles.size());
    }

private:
    // visits each owned cell of every tile (tile's own per-cell bucket + AMR cell),
    // reusing Tile::on_reachable_cells; LevelGhost skips non-clamp-owned cells
    template<ParticleType ptype, typename Fn>
    static void visit(ParticleArray_ const& particles, Fn&& fn)
    {
        auto const& tile_set = particles();
        for (auto const& tile : tile_set)
        {
            auto const& pc = tile();
            tile.template on_reachable_cells<ptype>([&](auto const& amr) {
                if constexpr (ptype == ParticleType::LevelGhost)
                    if (tile_set.at(amr) != &tile)
                        return;
                fn(pc, amr);
            });
        }
    }
};


template<typename ParticleArray_>
class ParticlePacker
{
    auto constexpr static dim = ParticleArray_::dimension;
    using ParticleArray_t     = ParticleArray_;
    using Particle_t          = ParticleDefaults<dim>::Particle_t;
    using Backend             = PackerBackend<ParticleArray_t::layout_mode, ParticleArray_t>;

    constexpr static Particle_t default_particle{};

public:
    static constexpr std::size_t n_keys = 5;

    ParticlePacker(ParticleArray_t const& particles)
        : particles_{particles}
    {
    }


    NO_DISCARD static auto get(Particle_t const& particle)
    {
        return std::forward_as_tuple(particle.weight_, particle.charge_, particle.iCell_,
                                     particle.delta_, particle.v_);
    }

    static constexpr auto empty() { return get(default_particle); }


    // sometimes we use this to infer the size of an ParticleArray
    // could be "charge" either
    NO_DISCARD static auto arbitrarySingleValueKey() { return "weight"; }

    NO_DISCARD static auto keys() { return packer_keys(); }

    NO_DISCARD auto get(std::size_t i) const { return get(particles_[i]); }
    NO_DISCARD bool hasNext() const { return it_ < particles_.size(); }
    NO_DISCARD auto next() { return get(it_++); }

    template<ParticleType ptype = ParticleType::Domain>
    NO_DISCARD static std::size_t size(ParticleArray_t const& particles)
    {
        return Backend::template size<ptype>(particles);
    }

    template<ParticleType ptype = ParticleType::Domain, typename SoAParticles_t>
    void pack(SoAParticles_t& copy)
    {
        Backend::template pack<ptype>(particles_, copy);
    }

private:
    ParticleArray_t const& particles_;
    std::size_t it_ = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLE_PACKER_H */
