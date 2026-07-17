#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include "core/def.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"
#include "core/data/particles/particle_array_sorter.hpp"
#include "core/data/particles/particle_array_selector.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"

#include "core/utilities/equality.hpp"
#include "core/utilities/monitoring.hpp"
#include "core/utilities/span.hpp"

#include <tuple>
#include <utility>
#include <sstream>
#include <type_traits>

namespace PHARE::core
{

template<auto opts /* defaulted in details header */>
class ParticleArray : public ResolvedParticleArray_t<opts>
{
    using This      = ParticleArray<opts>;
    using internals = ParticleArrayResolver<opts>;


    template<typename... Args>
    bool consteval static require()
    {
        using Tup = std::tuple<Args...>;

        bool constexpr base = std::is_constructible_v<Super, Args&&...>;
        if constexpr (std::tuple_size_v<Tup> > 0)
        {
            bool constexpr isself
                = std::is_same_v<std::decay_t<std::tuple_element_t<0, Tup>>, This>;
            // static_assert(!isself);
            return !isself and base;
        }
        else
            return base;
    }

public:
    using Super      = ResolvedParticleArray_t<opts>;
    using value_type = ParticleDefaults<opts.dim>::Particle_t;
    using view_t     = ParticleArray<opts.with_storage(StorageMode::SPAN)>;

    auto static constexpr options      = opts;
    auto static constexpr dimension    = opts.dim;
    auto static constexpr alloc_mode   = opts.alloc_mode;
    auto static constexpr layout_mode  = opts.layout_mode;
    auto static constexpr storage_mode = opts.storage_mode;
    auto static constexpr type_id      = internals::type_id;
    // auto static constexpr is_mapped    = internals::is_mapped; // torm
    auto static inline mon = MemoryMonitor{std::string{type_id}};

    std::string static id() { return std::string{type_id}; }

    ParticleArray(ParticleArray&& that)
        : Super{std::forward<Super>(that)}
    {
        mon.move();
    }
    ParticleArray(ParticleArray const& that)
        : Super{that}
    {
        mon.copy();
    }

    ParticleArray& operator=(ParticleArray&& that)
    {
        mon.move_assign();
        super() = std::move(that.super());
        return *this;
    }
    ParticleArray& operator=(ParticleArray const& that)
    {
        mon.copy_assign();
        super() = that.super();
        return *this;
    }

    template<typename... Args>
    ParticleArray(Args&&... args)
        requires(require<Args...>())
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
        mon.create();
    }



    auto view() { return view_t{*this}; }
    auto view() const { return view_t{*this}; }

    auto view(std::size_t i) // to take only i particles and ignore the rest
    {
        view_t v{*this};
        v.super().resize(i);
        return v;
    }

    auto view(std::size_t const start, std::size_t const size)
    {
        return view_t{*this, start, size};
    }

    auto operator*() { return view(); }
    auto operator*() const { return view(); }

    Super& super() _PHARE_ALL_FN_ { return *this; }
    Super const& super() const { return *this; }

    template<auto _opts>
    friend std::ostream& operator<<(std::ostream& out, ParticleArray<_opts> const&);
};

template<std::size_t dim>
using AoSParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>;


template<std::size_t dim>
using AoSMappedParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoSMapped}>;

template<std::size_t dim>
using SoAParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::SoA}>;


template<auto opts>
std::ostream& operator<<(std::ostream& out, ParticleArray<opts> const& arr)
{
    if constexpr (ParticleArray<opts>::layout_mode == LayoutMode::SoAPC)
    {
    }
    else
        for (auto const& p : arr)
            out << p.copy();
    return out;
}




template<auto opts>
void empty(ParticleArray<opts>& array)
{
    array.clear();
}


template<auto opts>
void swap(ParticleArray<opts>& array1, ParticleArray<opts>& array2)
{
    array1.swap(array2);
}

template<typename P0, typename P1>
EqualityReport particle_compare(P0 const& p0, P1 const& p1, std::size_t const i = 0,
                                double const atol = 1e-15)
{
    std::string idx = std::to_string(i);
    if (p0.iCell() != p1.iCell())
        return EqualityReport{false, "icell mismatch at index: " + idx, i};

    if (!float_equals(p0.v(), p1.v(), atol))
        return EqualityReport{false, "v mismatch at index: " + idx, i};
    if (!float_equals(p0.delta(), p1.delta(), atol))
        return EqualityReport{false, "delta mismatch at index: " + idx, i};


    return EqualityReport{true};
}


template<typename P0, typename P1>
EqualityReport particles_equals(P0 const& ref, P1 const& cmp, double const atol = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ref.size()) + " vs "
                                         + std::to_string(cmp.size())};

    auto rit      = ref.begin();
    auto cit      = cmp.begin();
    std::size_t i = 0;

    for (; rit != ref.end(); ++rit, ++cit, ++i)
        if (auto const eq = particle_compare(*rit, *cit, i, atol); !eq)
            return eq;

    return EqualityReport{true};
}

template<auto o>
EqualityReport operator==(ParticleArray<o> const& p0, ParticleArray<o> const& p1)
{
    auto report = particles_equals(p0, p1);
    if (!report)
    {
        PHARE_LOG_LINE_STR(p0[report.idx].copy());
        PHARE_LOG_LINE_STR(p1[report.idx].copy());
    }
    return report;
}

template<auto o>
EqualityReport operator==(ParticleArray<o> const& p0, std::vector<Particle<o>> const& p1)
{
    return particles_equals(p0, p1);
}



template<typename ParticleArray_t>
auto constexpr base_layout_type()
{
    using enum LayoutMode;
    auto constexpr layout_mode = ParticleArray_t::layout_mode;
    if constexpr (any_in(layout_mode, AoS, AoSMapped, AoSPC, AoSTS, AoSCMTS, AoSPCTS))
        return AoS;
    return SoA;
}



template<auto o>
void per_particle(ParticleArray<o> const& particles, auto const fn)
{
    using enum LayoutMode;
    if constexpr (is_tiled(o.layout_mode))
        for (auto& tile : particles())
            for (auto& p : tile())
                fn(p);
    else
        for (auto& p : particles)
            fn(p);
}

template<auto o>
void check_particles(ParticleArray<o> const& particles, [[maybe_unused]] bool const print = false)
{
    using enum LayoutMode;

    PHARE_DEBUG_DO({
        if constexpr (is_tiled(o.layout_mode))
            particles.check();
    })
}


template<auto o>
void check_particles_views(ParticleArray<o>& particles, [[maybe_unused]] bool const print = false)
{
    using enum LayoutMode;

    PHARE_DEBUG_DO({ check_particles(*particles); })
}

template<auto o>
void particle_array_domain_is_valid(ParticleArray<o> const& particles, auto const& domain_box)
{
    std::size_t in_domain_box = 0, not_in_domain_box = 0;

    if constexpr (o.layout_mode == LayoutMode::AoSPCTS)
    {
        for (auto const& tile : particles())
        {
            auto const& pc = tile();

            // every particle bucketed anywhere in this tile's ghost box (domain +
            // halo) must actually sit inside the tile itself, and must be filed
            // under the same cell its own iCell() maps to
            for (auto const& bix : pc.ghost_box())
            {
                auto const local           = pc.local_cell(bix);
                auto const& cell_particles = pc(local);
                for (std::size_t i = 0; i < cell_particles.size(); ++i)
                {
                    auto const& p = cell_particles[i];
                    if (not isIn(p, tile))
                    {
                        std::ostringstream oss;
                        oss << "particle_array_domain is not valid: particle outside tile"
                            << " storage_mode=" << static_cast<int>(o.storage_mode)
                            << " tile_box=" << tile << " bix=" << Point{bix}
                            << " local=" << Point{local} << " iCell=" << Point{p.iCell()};
                        throw std::runtime_error(oss.str());
                    }
                    if (not array_equals(pc.local_cell(p.iCell()), local))
                    {
                        std::ostringstream oss;
                        oss << "particle_array_domain is not valid: particle iCell does not "
                               "match its per-cell bucket"
                            << " storage_mode=" << static_cast<int>(o.storage_mode)
                            << " tile_box=" << tile << " bix=" << Point{bix}
                            << " local=" << Point{local} << " iCell=" << Point{p.iCell()};
                        throw std::runtime_error(oss.str());
                    }
                }
            }
        }
    }
    else if constexpr (is_tiled(o.layout_mode))
    {
        for (auto const& tile : particles())
            for (auto const& p : tile())
                if (not isIn(p, tile))
                    throw std::runtime_error("particle_array_domain is not valid");
    }

    per_particle(particles, [&](auto const& p) {
        if (isIn(p, domain_box))
            ++in_domain_box;
        else
            ++not_in_domain_box;
    });

    if (not(not_in_domain_box == 0 and in_domain_box == particles.size()))
        throw std::runtime_error("Invalid particles");

    // recurse once into the SPAN view so both storage sides are checked
    if constexpr (o.storage_mode == StorageMode::VECTOR)
        particle_array_domain_is_valid(*particles, domain_box);
}

template<auto o>
void particle_array_ghost_is_valid(ParticleArray<o> const& particles, auto const& domain_box,
                                   auto const& ghost_box)
{
    std::size_t in_ghost_layer = 0, not_in_ghost_layer = 0, outside_gb = 0;

    per_particle(particles, [&](auto const& p) {
        auto const ingb = isIn(p, ghost_box);
        auto const indb = isIn(p, domain_box);

        if (!ingb)
            ++outside_gb;
        else if (ingb and not indb)
            ++in_ghost_layer;
        else
            ++not_in_ghost_layer;
    });

    if constexpr (o.layout_mode == LayoutMode::AoSPCTS)
    {
        for (auto const& tile : particles())
        {
            auto const& pc = tile();

            // every particle bucketed anywhere in this tile's ghost box (domain +
            // halo) must actually sit in the halo, i.e. outside the patch's own
            // domain box, and must be filed under the same cell its own iCell() maps to
            for (auto const& bix : pc.ghost_box())
            {
                auto const local           = pc.local_cell(bix);
                auto const& cell_particles = pc(local);
                for (std::size_t i = 0; i < cell_particles.size(); ++i)
                {
                    auto const& p = cell_particles[i];
                    if (isIn(p, domain_box))
                    {
                        std::ostringstream oss;
                        oss << "particle_array_ghost is not valid: particle inside domain box"
                            << " storage_mode=" << static_cast<int>(o.storage_mode)
                            << " tile_box=" << tile << " bix=" << Point{bix}
                            << " local=" << Point{local} << " iCell=" << Point{p.iCell()};
                        throw std::runtime_error(oss.str());
                    }
                    if (not array_equals(pc.local_cell(p.iCell()), local))
                    {
                        std::ostringstream oss;
                        oss << "particle_array_ghost is not valid: particle iCell does not "
                               "match its per-cell bucket"
                            << " storage_mode=" << static_cast<int>(o.storage_mode)
                            << " tile_box=" << tile << " bix=" << Point{bix}
                            << " local=" << Point{local} << " iCell=" << Point{p.iCell()};
                        throw std::runtime_error(oss.str());
                    }
                }
            }
        }
    }
    else if constexpr (is_tiled(o.layout_mode))
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto const& tile = particles()[tidx];
            for (std::size_t i = 0; i < tile().size(); ++i)
                if (isIn(tile()[i], domain_box))
                    throw std::runtime_error("particle_array_ghost is not valid");
        }
    }

    if (not(outside_gb == 0 and not_in_ghost_layer == 0 and in_ghost_layer == particles.size()))
        throw std::runtime_error("Invalid particles");
}


// resolves to the final, most nested ParticleArray-like type for a given layout:
// itself if flat, a tile's particles if tiled, a cell's particles if per-cell
// (recursing once more into the tile's own per-cell particles for AoSPCTS).
template<auto o>
auto constexpr chunk_type_helper()
{
    using enum LayoutMode;
    if constexpr (any_in(o.layout_mode, AoSPC, SoAPC))
        return std::type_identity<typename ParticleArray<o>::per_cell_particles>{};
    else if constexpr (o.layout_mode == AoSPCTS)
    {
        using PerTile = typename ParticleArray<o>::per_tile_particles;
        using CellIt  = decltype(std::declval<PerTile&>().local_box().begin());
        return std::type_identity<std::remove_reference_t<decltype(std::declval<PerTile&>()(
            *std::declval<CellIt&>()))>>{};
    }
    else if constexpr (any_in(o.layout_mode, AoSTS, AoSCMTS, SoATS, SoAVXTS))
        return std::type_identity<typename ParticleArray<o>::per_tile_particles>{};
    else
        return std::type_identity<ParticleArray<o>>{};
}

template<auto o>
using ParticleArrayChunk_t = typename decltype(chunk_type_helper<o>())::type;


template<typename ParticleArray_t>
auto constexpr final_nested_type_helper()
{
    return chunk_type_helper<ParticleArray_t::options>();
}

template<typename ParticleArray_t>
using MostNestedParticleArray_t =
    typename decltype(final_nested_type_helper<ParticleArray_t>())::type;


// a view over a contiguous run of tiles that dereferences straight through to
// each tile's particles, so it can be range-for'd like a Span<Chunk> without
// copying anything out of the tiles themselves
template<typename Tile_t, typename Chunk_t>
class TileChunkSpan
{
public:
    class iterator
    {
    public:
        iterator(Tile_t* ptr, Tile_t* end)
            : ptr_{ptr}
            , end_{end}
        {
        }

        Chunk_t& operator*() const
        {
            assert(ptr_ != end_);
            return (*ptr_)();
        }

        iterator& operator++()
        {
            assert(ptr_ != end_); // don't walk the tile pointer past its own array
            ++ptr_;
            return *this;
        }

        bool operator!=(iterator const& that) const { return ptr_ != that.ptr_; }

    private:
        Tile_t* ptr_;
        Tile_t* end_;
    };

    TileChunkSpan(Tile_t* tiles, std::size_t size)
        : tiles_{tiles}
        , size_{size}
    {
    }

    auto begin() const { return iterator{tiles_, tiles_ + size_}; }
    auto end() const { return iterator{tiles_ + size_, tiles_ + size_}; }

private:
    Tile_t* tiles_;
    std::size_t size_;
};


// like TileChunkSpan, but keeps the tile itself alongside its particles, for
// callers that also need tile metadata (box, layout):
//   for (auto [tile, particles] : enumerate_tiles(ps)) { ... }
template<typename Tile_t, typename Chunk_t>
class TileSpan
{
public:
    class iterator
    {
    public:
        iterator(Tile_t* ptr, Tile_t* end)
            : ptr_{ptr}
            , end_{end}
        {
        }

        std::pair<Tile_t&, Chunk_t&> operator*() const
        {
            assert(ptr_ != end_);
            return {*ptr_, (*ptr_)()};
        }

        iterator& operator++()
        {
            assert(ptr_ != end_);
            ++ptr_;
            return *this;
        }

        bool operator!=(iterator const& that) const { return ptr_ != that.ptr_; }

    private:
        Tile_t* ptr_;
        Tile_t* end_;
    };

    TileSpan(Tile_t* tiles, std::size_t size)
        : tiles_{tiles}
        , size_{size}
    {
    }

    auto begin() const { return iterator{tiles_, tiles_ + size_}; }
    auto end() const { return iterator{tiles_ + size_, tiles_ + size_}; }

private:
    Tile_t* tiles_;
    std::size_t size_;
};


// enumerates at the tile level (one level of unwrapping only, unlike
// enumerate() which flattens all the way to leaf chunks) so callers that need
// tile metadata (box, layout) alongside its particles can get both:
//   for (auto [tile, particles] : enumerate_tiles(ps)) { tile.layout(); ... }
template<auto o>
auto enumerate_tiles(ParticleArray<o>& particles)
{
    using enum LayoutMode;
    static_assert(any_in(o.layout_mode, AoSTS, AoSCMTS, SoATS, SoAVXTS, AoSPCTS),
                  "enumerate_tiles() only makes sense for tiled layouts");

    using Chunk  = typename ParticleArray<o>::per_tile_particles;
    using Tile_t = std::remove_pointer_t<decltype(particles().data())>;
    return TileSpan<Tile_t, Chunk>{particles().data(), particles().size()};
}


// for AoSPCTS: a tile's own particles are themselves per-cell, so this
// iterates a TileChunkSpan of tiles (yielding each tile's per-cell container)
// and, for each one, walks that container's own cells (domain + halo) --
// flattening both levels of nesting down to the leaf particle chunks, with no
// allocation.
template<typename Tile_t, typename PerTileContainer_t, typename Chunk_t>
class PerCellTileChunkSpan
{
    using Outer   = TileChunkSpan<Tile_t, PerTileContainer_t>;
    using OuterIt = typename Outer::iterator;
    using CellBox = decltype(std::declval<PerTileContainer_t&>().local_box());
    using CellIt  = decltype(std::declval<CellBox&>().begin());

public:
    class iterator
    {
    public:
        iterator(OuterIt outer_it, OuterIt outer_end)
            : outer_it_{outer_it}
            , outer_end_{outer_end}
        {
            enter_tile();
        }

        Chunk_t& operator*() const
        {
            assert(outer_it_ != outer_end_);
            return (*outer_it_)(*cell_it_);
        }

        iterator& operator++()
        {
            assert(outer_it_ != outer_end_);
            ++cell_it_;
            if (*cell_it_ == *cell_end_)
            {
                ++outer_it_;
                enter_tile();
            }
            return *this;
        }

        bool operator!=(iterator const& that) const { return outer_it_ != that.outer_it_; }

    private:
        // advance outer_it_ to the next tile with at least one cell (full local box:
        // domain + halo, so ghost-content arrays are enumerated too, not just domain
        // ones), positioning cell_it_/cell_end_ on it; stops at outer_end_ if none left
        void enter_tile()
        {
            for (; outer_it_ != outer_end_; ++outer_it_)
            {
                auto& pc  = *outer_it_;
                cell_box_ = pc.local_box();
                cell_it_  = cell_box_.begin();
                cell_end_ = cell_box_.end();
                if (cell_it_ != cell_end_)
                    return;
            }
        }

        OuterIt outer_it_, outer_end_;
        CellBox cell_box_{};
        // box_iterator has no default ctor; seed both from cell_box_ (declared just
        // above) -- enter_tile() overwrites all three with real values regardless
        CellIt cell_it_{cell_box_.begin()};
        CellIt cell_end_{cell_box_.end()};
    };

    PerCellTileChunkSpan(Tile_t* tiles, std::size_t size)
        : outer_{tiles, size}
    {
    }

    auto begin() const { return iterator{outer_.begin(), outer_.end()}; }
    auto end() const { return iterator{outer_.end(), outer_.end()}; }

private:
    Outer outer_;
};


// enumerates the leaf-level particle chunks of a (possibly tiled/per-cell)
// ParticleArray, so callers don't need to know if the array is tiled to loop
// over its particles. Never allocates: it's a thin view over storage the
// particle array already owns (the tile/per-cell array itself, or the whole
// array when it's already flat).
template<auto o>
auto enumerate(ParticleArray<o>& particles)
{
    using enum LayoutMode;
    using Chunk = ParticleArrayChunk_t<o>;

    if constexpr (any_in(o.layout_mode, AoSPC, SoAPC))
        return Span<Chunk>{particles.data(), particles.size()};
    else if constexpr (o.layout_mode == AoSPCTS)
    {
        using Tile_t  = std::remove_pointer_t<decltype(particles().data())>;
        using PerTile = typename ParticleArray<o>::per_tile_particles;
        return PerCellTileChunkSpan<Tile_t, PerTile, Chunk>{particles().data(), particles().size()};
    }
    else if constexpr (any_in(o.layout_mode, AoSTS, AoSCMTS, SoATS, SoAVXTS))
    {
        using Tile_t = std::remove_pointer_t<decltype(particles().data())>;
        return TileChunkSpan<Tile_t, Chunk>{particles().data(), particles().size()};
    }
    else // flat: AoS, AoSMapped, SoA, SoAVX
        return Span<Chunk>{&particles, 1};
}


} // namespace PHARE::core


#endif
