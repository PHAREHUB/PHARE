#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CELL_MAPPED_TILE_SET_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CELL_MAPPED_TILE_SET_HPP

#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

// one flat, cell-mapped (AoSMapped) particle vector per tile: within-tile moves only touch
// the tile's cellmap, only cross-tile moves copy/delete particles
template<typename Particles_t>
class MappedParticlesTile : public Box<std::int32_t, Particles_t::dimension>
{
    static constexpr auto dim = Particles_t::dimension;
    using Super               = Box<std::int32_t, dim>;
    using This                = MappedParticlesTile<Particles_t>;

public:
    static constexpr auto dimension = dim;

    // the tile's cellmap spans its ghost box: level ghost / off-tile particles are mapped too
    MappedParticlesTile(Super const& box, std::size_t const ghost_cells)
        : Super{box}
        , particles{grow(box, ghost_cells)}
    {
    }

    MappedParticlesTile(MappedParticlesTile const&)            = default;
    MappedParticlesTile(MappedParticlesTile&&)                 = default;
    MappedParticlesTile& operator=(MappedParticlesTile const&) = default;
    MappedParticlesTile& operator=(MappedParticlesTile&&)      = default;

    template<typename OtherParticles_t>
    MappedParticlesTile(MappedParticlesTile<OtherParticles_t>& other)
        : Super{other}
        , particles{other.particles}
    {
    }

    auto& operator()() { return particles; }
    auto& operator()() const { return particles; }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }

    template<ParticleType type>
    void on_reachable_cells(auto&& fn) const
    {
        Super const& tbox = *this;
        for (auto const& amr : tbox)
            fn(amr);
        if constexpr (type == ParticleType::LevelGhost)
            for (auto const& amr : particles.box()) // box is AoSMapped ghost box!
                if (!isIn(amr, tbox))
                    fn(amr);
    }

    Particles_t particles;

    // indices of particles registered by move_check as leaving for another tile's domain
    std::vector<std::size_t> leavers;
};


template<std::size_t dim>
class MappedTileSetBoxes
{
public:
    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }
    auto ghost_cells() const { return ghost_cells_; }

    auto local_cell(std::array<int, dim> const& icell) const
    {
        return as_local_cell(ghost_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const { return local_cell(icell.toArray()); }

    auto local_box() const
    {
        return box_from_zero_to_upper_minus_one(
            ghost_box_.shape().template toArray<std::uint32_t>());
    }
    auto local_box(Box<int, dim> const& from) const
    {
        return Box<std::uint32_t, dim>{local_cell(from.lower), local_cell(from.upper)};
    }

    // iterates [amr, local] cell pairs of box, as GridLayout::amr_lcl_idx
    auto amr_lcl_idx(Box<int, dim> const& box) const { return boxes_iterator{box, local_box(box)}; }

protected:
    MappedTileSetBoxes(Box<int, dim> const& box, std::size_t const ghost_cells)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
    {
    }

    std::size_t ghost_cells_;
    Box<int, dim> box_, ghost_box_;
};


template<typename Particles>
class MappedTileSetSpan : public MappedTileSetBoxes<Particles::dimension>
{
    using Base = MappedTileSetBoxes<Particles::dimension>;

public:
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;

    using This               = MappedTileSetSpan<Particles>;
    using Particle_t         = typename ParticleDefaults<dim>::Particle_t;
    using per_tile_particles = Particles;
    using Tile_t             = MappedParticlesTile<Particles>;

    template<typename MappedTileSetArray>
    MappedTileSetSpan(MappedTileSetArray& arr)
        : Base{arr.box(), arr.ghost_cells()}
        , size_{arr.size()}
        , particles_{resolve(arr)}
    {
    }

    auto size() const { return size_; }

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }

protected:
    template<typename MappedTileSetArray>
    auto static resolve(MappedTileSetArray& arr)
    {
        if constexpr (MappedTileSetArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
            return arr.particles_views_.make_view();
    }

    template<typename>
    friend class MappedTileSetSpan;

    std::size_t size_;
    TileSetView<Tile_t> particles_;

}; // MappedTileSetSpan


template<typename Particles>
class MappedTileSetVector : public MappedTileSetBoxes<Particles::dimension>
{
    using This = MappedTileSetVector<Particles>;
    using Base = MappedTileSetBoxes<Particles::dimension>;

    template<typename>
    friend class MappedTileSetSpan;

protected:
    using Base::box_;
    using Base::ghost_cells_;

public:
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using Particle_t                   = ParticleDefaults<dim>::Particle_t;
    using value_type                   = Particle_t;
    using PSpan_t                      = typename Particles::view_t;
    using per_tile_particles           = Particles;
    using Tile_t                       = MappedParticlesTile<Particles>;
    using SpnTile                      = MappedParticlesTile<PSpan_t>;

    MappedTileSetVector(box_t const& box, auto const ghost_cells)
        : Base{box, static_cast<std::size_t>(ghost_cells)}
    {
    }

    // particles_views_ is rebuilt from the new particles_ by its default member initializer
    MappedTileSetVector(MappedTileSetVector&& that)
        : Base{that}
        , particles_{std::move(that.particles_)}
        , total_size{that.total_size}
    {
    }

    MappedTileSetVector(MappedTileSetVector const& that)
        : Base{that}
        , particles_{that.particles_}
        , total_size{that.total_size}
    {
    }

    // the defaulted memberwise copy/move leaves particles_views_ stale, reconstruct instead
    MappedTileSetVector& operator=(MappedTileSetVector&& that)
    {
        if (this == &that)
            return *this;
        this->~MappedTileSetVector();
        new (this) MappedTileSetVector(std::move(that));
        return *this;
    }
    MappedTileSetVector& operator=(MappedTileSetVector const& that)
    {
        if (this == &that)
            return *this;
        this->~MappedTileSetVector();
        new (this) MappedTileSetVector(that);
        return *this;
    }

    auto size() const { return total_size; }

    void emplace_back(Particle_t const& p)
    {
        auto* tile = particles_.at(Point<int, dim>{p.iCell()});
        assert(tile);
        (*tile)().emplace_back(p); // maps it into the tile's cellmap
        ++total_size;
    }

    template<typename... Args>
    void emplace_back(double const weight, Args&&... args)
    {
        this->emplace_back(Particle_t{weight, args...});
    }

    void push_back(Particle_t const& p) { emplace_back(p); }

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }

    // the particle already carries its NEW cell, pt carries the OLD cell and the tile
    // physically holding it. Within-tile moves only update the tile's cellmap,
    // cross-tile moves are unmapped and registered as leavers for on_moved()
    template<auto particle_type>
    auto& move_check(auto const& pt, std::size_t const idx, auto const& particle);

    // applies move_check-registered cross-tile moves: each tile pulls its neighbours'
    // leavers it owns, then swap-pops its own leavers, then on_appended() to finalize
    template<auto type>
    void on_moved();

    // refreshes size/views after particles were added outside of emplace_back
    template<auto type = ParticleType::Domain>
    void on_appended();

    // removes every particle matching pred, handing each to on_erase first. Descending
    // swap-pop per tile, keeping the tile cellmap consistent without a full remap
    void erase_if(auto&& pred, auto&& on_erase);
    void erase_if(auto&& pred)
    {
        erase_if(pred, [](auto const&) {});
    }

    // removes the given mapped particles of tile ti (any order, no duplicates) by descending
    // swap-pop. No on_appended(): callers erase over several tiles, then call it once
    void erase_in_tile(std::size_t const ti, std::vector<std::size_t>& idxs);

    // the views' tile set never changes shape, only each tile's span needs refreshing
    void reset_views()
    {
        update_from([&](std::size_t const i) { return SpnTile{particles_[i]}; }, particles_views_);
    }

    void clear();

    auto& views() { return particles_views_; }
    auto& views() const { return particles_views_; }

    NO_DISCARD bool is_consistent() const;

protected:
    // tiles build from amrbox, but `.at()` function maps ghosts box
    GhostClampedTileSet<Tile_t, alloc_mode> particles_{box_, ghost_cells_, ghost_cells_};
    GhostClampedTileSet<SpnTile, alloc_mode> particles_views_
        = GhostClampedTileSet<SpnTile, alloc_mode>::make_from(
            [](auto& tile) -> auto& { return tile; }, particles_);

    // per tile, the other tiles a particle can reach it from in one step
    std::vector<std::vector<std::size_t>> neighbours_ = make_neighbours();

    std::size_t total_size = 0;

private:
    auto make_neighbours() const
    {
        std::vector<std::vector<std::size_t>> neighbours(particles_.size());
        for (std::size_t i = 0; i < particles_.size(); ++i)
            for (std::size_t j = 0; j < particles_.size(); ++j)
                if (i != j and grow(particles_[i], ghost_cells_) * particles_[j])
                    neighbours[i].push_back(j);
        return neighbours;
    }

}; // MappedTileSetVector<Particles>


template<typename Super_>
struct MappedTileSetParticles : public Super_
{
    using Super              = Super_;
    using This               = MappedTileSetParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_tile_particles = typename Super::per_tile_particles;

    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Super::size;

    MappedTileSetParticles(MappedTileSetParticles&&)                 = default;
    MappedTileSetParticles& operator=(MappedTileSetParticles&&)      = default;
    MappedTileSetParticles(MappedTileSetParticles const&)            = default;
    MappedTileSetParticles& operator=(MappedTileSetParticles const&) = default;

    template<typename... Args>
    MappedTileSetParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
        : Super{std::forward<Args>(args)...}
    {
    }

    // no begin()/end() here - tiled, use per_particle() instead

    auto nbr_particles_in(std::array<int, dimension> const arr) const
    {
        return (*this->particles_.at(Point<int, dimension>{arr}))().nbr_particles_in(arr);
    }

    auto nbr_particles_in(Box<int, dimension> const box) const
    {
        std::size_t n_particles = 0;
        for (auto const& cell : box)
            n_particles += nbr_particles_in(cell.toArray());
        return n_particles;
    }

    void print() const {}
    void check() const {}

}; // MappedTileSetParticles<Super>


template<typename Particles>
template<auto particle_type>
auto& MappedTileSetVector<Particles>::move_check(auto const& pt, std::size_t const idx,
                                                 auto const& particle)
{
    using enum ParticleType;
    static_assert(any_in(particle_type, Domain, LevelGhost));

    auto const& newcell = particle.iCell();
    if (array_equals(newcell, pt.icell))
        return *this; // old cell == new cell, no change required

    auto& tile = *particles_.at(pt.tile_cell);
    auto& ps   = tile();
    ps.map(pt.icell).remove(idx);

    bool const into_other_tile = isIn(newcell, box_) and not isIn(newcell, tile);
    bool leaves                = into_other_tile;
    if constexpr (particle_type == LevelGhost)
        if (not isIn(newcell, box_))
            leaves = not isIn(newcell, this->ghost_box())
                     or particles_.at(Point<int, dim>{newcell}) != &tile;

    if (leaves)
        tile.leavers.push_back(idx); // stays unmapped, removed in on_moved
    else
        ps.map(newcell).add(idx); // same tile, or patch leaver kept in the tile ghost box

    return *this;
}

template<typename Particles>
template<auto type>
void MappedTileSetVector<Particles>::on_moved()
{
    // pull: each tile is the only writer of its own vector and cellmap
    for (std::size_t ti = 0; ti < particles_.size(); ++ti)
    {
        auto& dst = particles_[ti];
        for (auto const ni : neighbours_[ti])
        {
            auto& src = particles_[ni];
            for (auto const idx : src.leavers)
            {
                auto const cell = Point<int, dim>{src()[idx].iCell()};
                if constexpr (type == ParticleType::LevelGhost)
                    if (not isIn(cell, this->ghost_box())) // left the level ghost box
                        continue;
                if (particles_.at(cell) == &dst)
                    dst().emplace_back(src()[idx]); // maps it
            }
        }
    }

    // descending swap-pop keeps pending indices valid: the last particle is never a
    // pending leaver, and is still mapped at its current cell, so remap it to its new index
    for (auto& tile : particles_)
    {
        auto& ps = tile();
        auto& ls = tile.leavers;
        std::sort(ls.begin(), ls.end(), std::greater<>());
        for (auto const idx : ls)
        {
            auto const last = ps.size() - 1;
            if (idx != last)
            {
                ps.map(ps[last].iCell()).updateIndex(last, idx);
                ps.assign(last, idx);
            }
            ps.pop_back();
        }
        ls.clear();
    }

    on_appended<type>();
}

template<typename Particles>
void MappedTileSetVector<Particles>::erase_if(auto&& pred, auto&& on_erase)
{
    for (auto& tile : particles_)
    {
        auto& ps = tile();
        for (std::size_t i = ps.size(); i-- > 0;)
        {
            if (!pred(ps[i]))
                continue;
            on_erase(ps[i]);
            ps.map(ps[i].iCell()).remove(i);
            auto const last = ps.size() - 1;
            if (i != last) // last was already checked and kept
            {
                ps.map(ps[last].iCell()).updateIndex(last, i);
                ps.assign(last, i);
            }
            ps.pop_back();
        }
    }
    on_appended();
}

template<typename Particles>
void MappedTileSetVector<Particles>::erase_in_tile(std::size_t const ti,
                                                   std::vector<std::size_t>& idxs)
{
    // descending keeps pending indices valid: the last particle is never a pending one
    auto& ps = particles_[ti]();
    std::sort(idxs.begin(), idxs.end(), std::greater<>());
    for (auto const i : idxs)
    {
        ps.map(ps[i].iCell()).remove(i);
        auto const last = ps.size() - 1;
        if (i != last)
        {
            ps.map(ps[last].iCell()).updateIndex(last, i);
            ps.assign(last, i);
        }
        ps.pop_back();
    }
}

template<typename Particles>
template<auto type>
void MappedTileSetVector<Particles>::on_appended()
{
    total_size = 0;
    for (auto const& tile : particles_)
        total_size += tile().size();

    reset_views();
}

template<typename Particles>
void MappedTileSetVector<Particles>::clear()
{
    for (auto& tile : particles_)
        tile().clear(); // clears the tile's cellmap too
    total_size = 0;

    reset_views();
}

template<typename Particles>
NO_DISCARD bool MappedTileSetVector<Particles>::is_consistent() const
{
    std::size_t n = 0;
    for (auto const& tile : particles_)
    {
        if (!tile().is_consistent())
            return false;
        n += tile().size();
    }
    return n == total_size;
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CELL_MAPPED_TILE_SET_HPP */
