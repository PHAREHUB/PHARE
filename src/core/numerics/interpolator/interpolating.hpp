#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "core/data/field/field_tiles.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "interpolator.hpp"

namespace PHARE::core
{

// visits, for tile tidx, the particles of every cell its level ghost halo reaches, read
// from that cell's clamp-owner tile (TileSet::tag_cells_). Level ghosts are never in the
// patch domain, so only halo cells outside of it are visited: none for interior tiles
template<typename Particles>
void on_reachable_level_ghosts(Particles const& particles, std::size_t const tidx, auto&& fn)
{
    auto const& patch_box    = particles.box();
    auto const ghost_cells   = patch_box.lower[0] - particles.ghost_box().lower[0];
    auto const& tile         = particles()[tidx];
    auto const tile_halo     = grow(*tile, ghost_cells);
    auto const halo_in_patch = tile_halo * patch_box;
    if (halo_in_patch and *halo_in_patch == tile_halo)
        return; // interior tile

    for (auto const& amr : tile_halo)
    {
        if (isIn(amr, patch_box))
            continue;
        auto const& owner = (*particles().at(amr))();
        if constexpr (Particles::layout_mode == LayoutMode::AoSCMTS)
            for (auto const idx : owner.map(amr))
                fn(makeRange(owner, idx, idx + 1));
        else
            fn(owner(owner.local_cell(amr)));
    }
}

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false,
         typename Interpolator_t = Interpolator<dim, interpOrder, atomic_ops>>
class Interpolating
{
public:
    template<typename Particles>
    inline void operator()(Particles const& particles, auto& rhoP, auto& flux, auto const& layout,
                           double coef = 1.)
    {
        particleToMesh(particles, layout, rhoP, flux, coef);
    }

    // LevelGhost: level ghosts live only in their clamp-owner tile, and their deposits
    // happen after the tile reduction, so each tile gathers (read-only) every level ghost
    // its halo reaches from the owning tile
    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(is_tiled(Particles::layout_mode))
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
            particleToMesh_tile<type>(tidx, particles, rhoP, flux, coef);
    }

    // only writes the moments of tile tidx, distinct tiles can be deposited concurrently
    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh_tile(std::size_t const tidx, Particles const& particles, auto& rhoP,
                             auto& flux, double coef = 1.)
        requires(is_tiled(Particles::layout_mode))
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        auto [rhop, F]      = tiles_at(tidx, rhoP, flux);
        auto const& rho_lay = tile_layout(rhoP, tidx);
        auto const deposit  = [&](auto const& cell_particles) {
            for (auto const& p : cell_particles)
                interp_.particleToMesh(p, rhop, F, rho_lay, coef);
        };

        auto& tile_particles = particles()[tidx]();
        if constexpr (type == ParticleType::LevelGhost)
            on_reachable_level_ghosts(particles, tidx, deposit);
        else if constexpr (Particles::layout_mode == LayoutMode::AoSCMTS)
            deposit(tile_particles);
        else
            for (auto const& bix : tile_particles.local_box())
                deposit(tile_particles(bix));
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& flux,
                        double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSMapped)
    {
        interp_(particles, rhoP, flux, layout, coef);
    }

    Interpolator_t interp_;
};

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false>
struct MomentumTensorInterpolating
{
    using Interpolator_t = MomentumTensorInterpolator<dim, interpOrder, atomic_ops>;

    Interpolator_t interp_;

public:
    template<typename Particles_t>
    inline void operator()(Particles_t& particles, auto& momentumTensor, auto const& layout,
                           double mass = 1.)
    {
        particleToMesh(particles, momentumTensor, layout, mass);
    }

    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSMapped)
    {
        interp_(particles, momentumTensor, layout, mass);
    }

    // see Interpolating::particleToMesh
    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(is_tiled(Particles_t::layout_mode))
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt            = tile_at(momentumTensor, tidx);
            auto const& mt_lay = tile_layout(momentumTensor, tidx);
            auto const deposit
                = [&](auto const& cell_particles) { interp_(cell_particles, mt, mt_lay, mass); };

            auto& tile_particles = particles()[tidx]();
            if constexpr (type == ParticleType::LevelGhost)
                on_reachable_level_ghosts(particles, tidx, deposit);
            else if constexpr (Particles_t::layout_mode == LayoutMode::AoSCMTS)
                deposit(tile_particles);
            else
                for (auto const& bix : tile_particles.local_box())
                    deposit(tile_particles(bix));
        }
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
