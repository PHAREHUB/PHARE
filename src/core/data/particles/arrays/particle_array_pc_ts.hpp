#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PER_CELL_TILE_SET_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PER_CELL_TILE_SET_HPP

#include "core/def.hpp"
#include "core/operators.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_translation_tracking.hpp"

namespace PHARE::core
{

template<typename PerCellParticles>
class PCParticlesTile : public Box<std::int32_t, PerCellParticles::dimension>
{
    static constexpr auto dim = PerCellParticles::dimension;
    using Super               = Box<std::int32_t, dim>;
    using This                = PCParticlesTile<PerCellParticles>;

public:
    static constexpr auto dimension = dim;

    PCParticlesTile(Super const& box, std::size_t const ghost_cells)
        : Super{box}
        , particles{box, ghost_cells}
    {
    }

    PCParticlesTile(PCParticlesTile const&)            = default;
    PCParticlesTile(PCParticlesTile&&)                 = default;
    PCParticlesTile& operator=(PCParticlesTile const&) = default;
    PCParticlesTile& operator=(PCParticlesTile&&)      = default;

    template<typename OtherPerCell>
    PCParticlesTile(PCParticlesTile<OtherPerCell>& other)
        : Super{other}
        , particles{other.particles}
    {
    }

    auto& operator()() { return particles; }
    auto& operator()() const { return particles; }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }

    auto& link(std::size_t const idx) { return _links[idx]; }
    auto& links() { return _links; }
    auto& links() const { return _links; }

    // every cell this tile physically reaches: its domain box, plus - for LevelGhost -
    // its own ghost halo (unfiltered: adjacent tiles' halos legitimately overlap).
    // callers wanting a deduped subset filter fn themselves
    template<ParticleType type>
    void on_reachable_cells(auto&& fn) const
    {
        Super const& tbox = *this;
        for (auto const& amr : tbox)
            fn(amr);
        if constexpr (type == ParticleType::LevelGhost)
            for (auto const& amr : particles.ghost_box())
                if (!isIn(amr, tbox))
                    fn(amr);
    }

    PerCellParticles particles;

private:
    std::array<This*, 7> _links = ConstArray<This*, 7>(nullptr);
};


template<typename Particles>
struct PCCrossTileCopyDAO;


template<typename Particles>
class PCTileSetSpan : public ParticleTranslationTrackerSpan<Particles::dimension>
{
    using Base = ParticleTranslationTrackerSpan<Particles::dimension>;

    template<typename>
    friend struct PCCrossTileCopyDAO;

protected:
    using Base::add_into_;
    using Base::box_;
    using Base::cap_;
    using Base::gap_idx_;
    using Base::gaps_;
    using Base::ghost_box_;
    using Base::left_;
    using Base::size_;

public:
    auto static constexpr alloc_mode = Particles::alloc_mode;
    auto static constexpr dim        = Particles::dimension;

    using This               = PCTileSetSpan<Particles>;
    using lobox_t            = Box<std::uint32_t, dim>;
    using per_tile_particles = Particles;
    using Tile_t             = PCParticlesTile<Particles>;

    using Base::box;
    using Base::ghost_box;
    using Base::local_box;
    using Base::local_cell;
    using Base::size;

private:
    using locell_t = std::array<std::uint32_t, dim>;

    template<typename PCTileSetArray>
    auto resolve(PCTileSetArray& arr)
    {
        if constexpr (PCTileSetArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
            return arr.particles_views_.make_view();
    }

    template<typename PCTileSetArray>
    auto resolve_gaps(PCTileSetArray& arr)
    {
        if constexpr (PCTileSetArray::storage_mode == StorageMode::SPAN)
            return arr.gaps_;
        else
            return *arr.gap_views_;
    }

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;

    template<typename PCTileSetArray>
    PCTileSetSpan(PCTileSetArray& arr);

    auto size(std::size_t const& idx) const { return particles_.data()[idx]().size(); }
    auto size(locell_t const& icell) const { return cell_size_(icell); }

    auto& operator()() { return particles_; }
    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return (*particles_.at(cell))(); }
    auto& operator()(locell_t const& cell) const { return (*particles_.at(cell))(); }

    auto local_tile_cell(std::array<int, dim> const& cell) const
    {
        PHARE_ASSERT(particles_.at(local_cell(cell)));
        return local_cell((*particles_.at(local_cell(cell))).lower);
    }

    template<auto particle_type>
    auto& move_check(auto const& pt, std::size_t const idx, auto& particle);

    template<auto type, typename... Args>
    void sync(Args&&... args);

    void clear() {} // TODO

protected:
    template<auto type>
    void sync_tile_add_new(std::size_t const tidx);
    template<auto type>
    void sync_tile_rm_left(std::size_t const tidx);

    TileSetView<Tile_t> particles_;
    NdArrayView<dim, std::size_t> cell_size_;

}; // PCTileSetSpan


template<typename Particles>
class PCTileSetVector
    : public ParticleTranslationTracker<Particles::dimension, Particles::alloc_mode>
{
    using This = PCTileSetVector<Particles>;
    using Base = ParticleTranslationTracker<Particles::dimension, Particles::alloc_mode>;

    template<typename P>
    friend class PCTileSetSpan;

protected:
    using Base::add_into_;
    using Base::box_;
    using Base::cap_;
    using Base::cell_size_;
    using Base::gap_idx_;
    using Base::gap_views_;
    using Base::gaps_;
    using Base::ghost_box_;
    using Base::ghost_cells_;
    using Base::left_;
    using Base::total_size;

public:
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr layout_mode  = Particles::layout_mode;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = ParticleDefaults<dim>::Particle_t;
    using value_type                   = Particle_t;
    using PSpan_t                      = typename Particles::view_t;
    using per_tile_particles           = Particles;
    using Tile_t                       = PCParticlesTile<Particles>;
    using SpnTile                      = PCParticlesTile<PSpan_t>;
    using size_t_vector                = std::vector<std::size_t>;
    using Base::box;
    using Base::ghost_box;
    using Base::local_box;
    using Base::local_cell;
    using Base::zero_bookkeeping;

    PCTileSetVector(box_t const& box, auto const ghost_cells)
        : Base{box, ghost_cells}
    {
        TileSet<Tile_t, alloc_mode>::build_links(particles_);
        TileSet<SpnTile, alloc_mode>::build_links(particles_views_);
    }

    PCTileSetVector(PCTileSetVector&& that)
        : Base{that.box_, that.ghost_cells_}
        , particles_{std::move(that.particles_)}
    {
        on_appended(); // without std::swap does not work well - mirrors TileSetVector
    }

    PCTileSetVector(PCTileSetVector const& that)
        : Base{that.box_, that.ghost_cells_}
        , particles_{that.particles_.copy(TileSetter<dim>{that.box_, that.ghost_cells_})}
    {
        on_appended();
    }

    // not = default: see TileSetVector::operator= (particle_array_ts.hpp) -- the
    // defaulted memberwise copy/move leaves particles_views_ stale, so reconstruct
    // in place via the copy/move constructors instead, which call on_appended() themselves
    PCTileSetVector& operator=(PCTileSetVector&& that)
    {
        if (this == &that)
            return *this;
        this->~PCTileSetVector();
        new (this) PCTileSetVector(std::move(that));
        return *this;
    }
    PCTileSetVector& operator=(PCTileSetVector const& that)
    {
        if (this == &that)
            return *this;
        this->~PCTileSetVector();
        new (this) PCTileSetVector(that);
        return *this;
    }

    auto size() const { return total_size; }
    auto size(locell_t const& icell) const { return cell_size_(icell); }
    auto size(std::size_t const& idx) const { return cell_size_.data()[idx]; }

    void emplace_back(Particle_t const& p)
    {
        auto* tile = particles_.at(Point<int, dim>{p.iCell()});
        assert(tile);
        (*tile)().emplace_back(p);
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

    // refreshes recount/gap-sizing/views after particles were added without going through
    // move_check (a raw append) - see on_moved() for the post-move_check pipeline.
    template<auto type = ParticleType::Domain>
    void on_appended();

    template<auto type>
    void sync_moved(); // realloc + resize ahead of the view-side copies

    // applies move_check-registered moves: realloc/resize, then the span-side within-tile
    // adds + cross-tile copies + rm, then finalize (recount, size wall gaps, reset views).
    template<auto type>
    void on_moved(auto&&... args);

    void reset_views()
    {
        update_from([&](std::size_t const i) { return SpnTile{particles_[i]}; }, particles_views_);
        update_from([&](std::size_t const i) { return make_span(*(gaps_.data() + i)); },
                    gap_views_);
    }

    void clear();

    auto& views() { return particles_views_; }
    auto& views() const { return particles_views_; }

    NO_DISCARD bool is_consistent() const;

protected:
    // visit the AMR cells of a tile that sit within ghost_cells_ of its wall — the only
    // cells a particle can enter or leave the tile from in one step
    void on_tile_wall_cells(auto const& tile, auto&& fn) const;

    // visit the tile's ghost-box cells it clamp-owns (TileSet::tag_cells_) — the cells
    // level ghost particles live in and register their cross-tile traffic against
    void on_tile_owned_ghost_cells(auto const& tile, auto&& fn) const
    {
        Box<int, dim> const& tbox = tile;
        for (auto const& amr : tile().ghost_box())
            if (!isIn(amr, tbox) and particles_.at(amr) == &tile)
                fn(amr);
    }

    // tiles build from amrbox, but `.at()` function maps ghosts box
    TileSet<Tile_t, alloc_mode> particles_{TileSetter<dim>{box_, ghost_cells_}, ghost_cells_};
    TileSet<SpnTile, alloc_mode> particles_views_ = TileSet<SpnTile, alloc_mode>::make_from(
        [](auto& tile) -> auto& { return tile; }, TileSetter<dim>{box_, ghost_cells_}, particles_);

}; // PCTileSetVector<Particles>


template<typename Super_>
struct PCTileSetParticles : public Super_
{
    using Super              = Super_;
    using This               = PCTileSetParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_tile_particles = typename Super::per_tile_particles;

    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Super::size;

    PCTileSetParticles(PCTileSetParticles&&)                 = default;
    PCTileSetParticles& operator=(PCTileSetParticles&&)      = default;
    PCTileSetParticles(PCTileSetParticles const&)            = default;
    PCTileSetParticles& operator=(PCTileSetParticles const&) = default;

    template<typename... Args>
    PCTileSetParticles(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
        : Super{std::forward<Args>(args)...}
    {
    }

    // no begin()/end() here - AoSPCTS is tiled, so top-level particle iteration is not
    // supported (see ParticleArray::begin()/end() in particle_array.hpp); use
    // enumerate()/per_particle() instead.

    auto data() const { return static_cast<Particle_t const*>(nullptr); } // TODO
    auto data() { return static_cast<Particle_t*>(nullptr); }             // TODO

    // move_check lives on PCTileSetSpan — no stub here or it shadows the span's version
    auto nbr_particles_in(std::array<int, dimension> const arr) const
    {
        auto const& tile = *this->particles_.at(this->local_cell(arr));
        return tile()(tile().local_cell(arr)).size();
    }

    auto nbr_particles_in(Box<int, dimension> const box) const
    {
        std::size_t n_particles = 0;
        for (auto const& cell : box)
            n_particles += nbr_particles_in(cell.toArray());
        return n_particles;
    }

    auto max_size() const { return std::size_t{0}; /* TODO */ }

    void print() const {}
    void check() const {}

    template<typename T>
    struct index_wrapper;
    auto operator[](std::size_t const& s) { return index_wrapper<This>{this, s}; }
    auto operator[](std::size_t const& s) const { return index_wrapper<This const>{this, s}; }

}; // PCTileSetParticles<Super>


template<typename Particles>
struct PCCrossTileCopyDAO
{
    auto static constexpr dim = Particles::dim;
    using Tile                = PCParticlesTile<typename Particles::per_tile_particles>;

    // no ownership filter: each tile's move_check bookkeeping is private to itself,
    // so a non-owned duplicate ghost cell still has its own pending moves to flush
    template<auto type>
    void on_owned_cells(auto&& fn)
    {
        tile.template on_reachable_cells<type>(fn);
    }

    template<auto type>
    void copy_in();

    // both leaver lists of a cell are consumed by one merged rm pass on the per-cell
    // container — see sync_rm_left for why they cannot be two separate passes
    template<auto type>
    void rm_left()
    {
        auto& pc = tile();

        on_owned_cells<type>([&](auto const& amr) {
            auto const cix = ps.local_cell(amr);
            pc.sync_rm_left(pc.local_cell(amr), ps.gaps_(cix).data(), ps.gap_idx_(cix),
                            ps.left_(cix));
        });
    }

    Particles& ps;
    std::size_t src_tile_idx;
    Tile& tile = ps()[src_tile_idx];
};


template<typename Particles>
template<typename PCTileSetArray>
PCTileSetSpan<Particles>::PCTileSetSpan(PCTileSetArray& arr)
    : Base{{.gaps            = resolve_gaps(arr),
            .gap_idx         = arr.gap_idx_,
            .add_into        = arr.add_into_,
            .cap             = arr.cap_,
            .left            = arr.left_,
            .size            = arr.size(),
            .box             = arr.box_,
            .ghost_box       = arr.ghost_box(),
            .local_ghost_box = arr.local_box()}}
    , particles_{resolve(arr)}
    , cell_size_{arr.cell_size_}
{
}

template<typename Particles>
template<auto particle_type>
auto& PCTileSetSpan<Particles>::move_check(auto const& pt, std::size_t const idx, auto& particle)
{
    using enum ParticleType;
    static_assert(any_in(particle_type, Domain, LevelGhost));

    auto const& newcell = particle.iCell();

    if (array_equals(newcell, pt.icell))
        return *this; // old cell == new cell, no change required

    bool constexpr static ATOMIC = true;
    using Op                     = Operators<std::size_t, ATOMIC>;

    // register the departure against the tile-set cell — mirrors
    // TileSetParticles::move_check: only register here, the actual cross-tile copy
    // (or removal) happens later during sync
    auto const leave = [&]() {
        auto const old_lcl_cell = local_cell(pt.icell);

        auto& gidx      = gap_idx_(old_lcl_cell);
        auto const nidx = Op{gidx}.increment_return_old();
        auto& gaps      = gaps_(old_lcl_cell);
        assert(nidx < gaps.size());
        gaps[nidx] = idx;
    };

    if constexpr (particle_type == Domain)
    {
        auto& old_tile = *particles_.at(pt.tile_cell);
        if (isIn(newcell, box()) and not isIn(newcell, old_tile))
        {
            leave();
            Op{add_into_(local_cell(newcell))}.increment_return_old();
            return *this;
        }
    }
    else // LevelGhost
    {
        auto& cur_tile = *particles_.at(pt.tile_cell);
        if (isIn(newcell, box()) and not isIn(newcell, cur_tile))
        {
            cur_tile().register_particle_ejection(pt, idx);
            return *this;
        }

        cur_tile().template move_check<particle_type>(pt, idx, particle);
        return *this;
    }

    // Domain: cell change register: still owned by (or ghosted into) the current
    // tile. move_check's own box()/old_tile check above covers cross-tile hand-off;
    // this fallthrough covers "outside the domain, assumed still in this tile's
    // ghost box" and "inside the domain, inside this tile".
    (*particles_.at(pt.tile_cell))().template move_check<particle_type>(pt, idx, particle);

    return *this;
}


template<typename Particles>
template<auto type>
void PCTileSetVector<Particles>::on_appended()
{
    total_size = 0;
    for (auto& tile : particles_)
    {
        tile().template on_appended<type>(); // per-tile recount + gap sizing + views
        total_size += tile().size();
    }

    // only cells within ghost_cells_ of a tile wall can register cross-tile
    // departures; size their gap vectors so move_check's operator[] has room.
    // level ghost particles live in (and leave from) clamp-owned ghost cells, so
    // those are sized too — type-agnostic: the initial post-fill sync may run as
    // Domain on a level ghost array
    for (auto& tile : particles_)
    {
        auto const size_cell = [&](auto const& amr_cell) {
            auto const c   = local_cell(amr_cell);
            auto const& cs = tile()(tile().local_cell(amr_cell)).size();
            cell_size_(c)  = cs;
            if (auto& gaps = gaps_(c); gaps.size() < cs)
                gaps.resize(cs);
        };
        on_tile_wall_cells(tile, size_cell);
        on_tile_owned_ghost_cells(tile, size_cell);
    }

    reset_views();
}

template<typename Particles>
template<auto type>
void PCTileSetVector<Particles>::sync_moved() // realloc + resize ahead of the view-side copies
{
    for (auto& tile : particles_)
    {
        Box<int, dim> const& tbox = tile;
        auto& pc                  = tile();

        // fold in cross-tile traffic per cell. Domain: gate by clamp ownership -
        // it's shared bookkeeping, so exactly one tile must fold it in. LevelGhost:
        // always fold in - bookkeeping there is per-tile-private, and an ownership
        // gate would wrongly skip cells this tile's own ghost_box() reaches but a
        // different tile clamp-owns.
        for (auto const& amr : pc.ghost_box())
        {
            std::size_t in = 0, out = 0;
            if constexpr (type == ParticleType::LevelGhost)
            {
                auto const cix = local_cell(amr);
                in             = add_into_(cix);
                out            = gap_idx_(cix);
                add_into_(cix) = 0;
            }
            else if (isIn(amr, tbox) or particles_.at(amr) == &tile)
            {
                auto const cix = local_cell(amr);
                in             = add_into_(cix);
                out            = gap_idx_(cix);
                add_into_(cix) = 0;
            }
            pc.template sync_moved<type>(pc.local_cell(amr), in, out);
        }
        pc.reset_views(); // per-cell spans start the copy step at the pre-move sizes
    }
    reset_views();
}

template<typename Particles>
template<auto type>
void PCTileSetVector<Particles>::on_moved(auto&&... args)
{
    sync_moved<type>();
    PCTileSetSpan<PSpan_t>{*this}.template sync<type>(args...);

#if PHARE_DEBUG
    auto& vws = views();
    for (std::size_t i = 0; i < particles_.size(); ++i)
        for (auto const& bix : particles_[i]().local_box())
        {
            auto const vsize = vws[i]()(bix).size();
            auto const psize = particles_[i]()(bix).size();
            if (vsize != psize)
                PHARE_LOG_LINE_SS("on_moved mismatch: tile " << i << " cell " << bix << " view="
                                                             << vsize << " vector=" << psize);
        }
    for (std::size_t i = 0; i < particles_.size(); ++i)
        for (auto const& bix : particles_[i]().local_box())
            assert(vws[i]()(bix).size() == particles_[i]()(bix).size());
#endif

    on_appended<type>();
}

template<typename Particles>
void PCTileSetVector<Particles>::clear()
{
    for (auto& tile : particles_)
        tile().clear(); // resets each tile's own per-cell bookkeeping too
    total_size = 0;

    // this level's own cross-tile move_check bookkeeping (registered gaps/counts
    // against tile-set cells) is separate from each tile's - reset it the same way
    // the constructor does, or clear() leaves it stale just like the per-cell case.
    zero_bookkeeping();

    reset_views();
}

template<typename Particles>
NO_DISCARD bool PCTileSetVector<Particles>::is_consistent() const
{
    std::size_t n = 0;
    for (auto const& tile : particles_)
    {
        auto const& pc = tile();
        for (auto const& cell : pc.local_box())
            for (auto const& p : pc(cell.toArray()))
                if (pc.local_cell(p.iCell()) != cell.toArray())
                    return false;
        n += pc.size();
    }
    return n == total_size;
}

template<typename Particles>
void PCTileSetVector<Particles>::on_tile_wall_cells(auto const& tile, auto&& fn) const
{
    Box<int, dim> const& tbox = tile;
    auto const w              = static_cast<int>(ghost_cells_);

    bool whole = false;
    for (std::size_t d = 0; d < dim; ++d)
        whole |= tbox.shape()[d] <= 2 * w;

    if (whole) // too small for an interior, every cell is a wall cell
    {
        for (auto const& bix : tbox)
            fn(bix);
        return;
    }
    for (auto const& b : tbox.remove(shrink(tbox, w)))
        for (auto const& bix : b)
            fn(bix);
}

template<typename Particles>
template<auto type>
void PCCrossTileCopyDAO<Particles>::copy_in()
{
    auto& pc = tile();

    // within-tile movers (including off-domain movers headed for the tile's own
    // ghost layer) — registered on the tile's own per-cell gap lists
    on_owned_cells<type>([&](auto const& amr) { pc.sync_add_new(pc.local_cell(amr)); });

    // cross-tile leavers — registered against the tile-set cell; each lands in
    // whichever tile owns its particle's new cell
    on_owned_cells<type>([&](auto const& amr) {
        auto const cix     = ps.local_cell(amr);
        auto const& n_gaps = ps.gap_idx_(cix);
        if (!n_gaps)
            return;
        {
            auto& gaps = ps.gaps_(cix);
            pc.sort(gaps.data(), gaps.data() + n_gaps);
        }
        auto& src        = pc(pc.local_cell(amr));
        auto& left       = ps.left_(cix);
        auto const& gaps = ps.gaps_(cix);
        for (std::size_t i = 0; i < n_gaps; ++i)
        {
            auto const& gidx    = gaps[n_gaps - (1 + i)];
            auto const& newcell = src.iCell(gidx);

            if constexpr (type == ParticleType::LevelGhost)
                if (not isIn(newcell, ps.ghost_box()))
                { // left the level ghost box: no destination, rm_left deletes it
                    ++left;
                    continue;
                }

            auto& dst_tile = *ps().at(ps.local_cell(newcell));
            auto& dst_pc   = dst_tile();
            if constexpr (type == ParticleType::Domain)
            { // level ghost dst cells are clamp-owned: outside the dst tile box
                PHARE_ASSERT(not isIn(src[gidx], tile));
                PHARE_ASSERT(isIn(src[gidx], dst_tile));
            }
            [[maybe_unused]] bool const ok
                = dst_pc.append_from(dst_pc.local_cell(newcell), src, gidx);
            PHARE_ASSERT(ok); // capacity is exact-reserved in sync_moved
            ++left;
        }
        PHARE_ASSERT(left == n_gaps);
    });
}


template<typename Particles>
template<auto type, typename... Args>
void PCTileSetSpan<Particles>::sync(Args&&... args)
{
    PHARE_LOG_SCOPE(3, "PCTileSetSpan::sync(stream)");

    for (std::size_t tidx = 0; tidx < particles_.size(); ++tidx)
        sync_tile_add_new<type>(tidx);

    for (std::size_t tidx = 0; tidx < particles_.size(); ++tidx)
        sync_tile_rm_left<type>(tidx);
}


template<typename Particles>
template<auto type>
void PCTileSetSpan<Particles>::sync_tile_add_new(std::size_t const tidx)
{
    PCCrossTileCopyDAO<std::decay_t<decltype(*this)>>{*this, tidx}.template copy_in<type>();
}


template<typename Particles>
template<auto type>
void PCTileSetSpan<Particles>::sync_tile_rm_left(std::size_t const tidx)
{
    PCCrossTileCopyDAO<std::decay_t<decltype(*this)>>{*this, tidx}.template rm_left<type>();
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PER_CELL_TILE_SET_HPP */
