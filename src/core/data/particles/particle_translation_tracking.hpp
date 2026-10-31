#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_TRANSLATION_TRACKING_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_TRANSLATION_TRACKING_HPP

#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

// Shared per-cell bookkeeping for the owning (VECTOR) side of PerCellVector / TileSetVector /
// PCTileSetVector: the gap/arrival/capacity counters used to register and apply particle moves,
// plus the box/local-cell plumbing they're indexed by. Identical across all three layouts.
//
// What is NOT here on purpose: sync/sync_moved/on_appended/move_check and friends. Those walk
// genuinely different index spaces per layout (per-cell, per-tile via a representative slot,
// per-tile-wall-cell layered on a nested per-cell level) and stay implemented separately.
template<std::size_t dim, auto alloc_mode>
class ParticleTranslationTracker
{
public:
    using box_t   = Box<int, dim>;
    using lobox_t = Box<std::uint32_t, dim>;

    template<typename T>
    using nd_array_t = NdArrayVector<dim, T, /*c_order=*/true, alloc_mode>;

    using size_t_vector = std::vector<std::size_t>;

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }

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
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

protected:
    ParticleTranslationTracker(box_t const& box, std::size_t const ghost_cells)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
    {
        zero_bookkeeping();
    }

    void zero_bookkeeping()
    {
        cell_size_.zero();
        gap_idx_.zero();
        add_into_.zero();
        left_.zero();
        cap_.zero();
    }

    std::size_t ghost_cells_;
    Box<int, dim> box_, ghost_box_;

    nd_array_t<size_t_vector> gaps_{local_box().shape()};
    nd_array_t<Span<std::size_t>> gap_views_{local_box().shape()};
    nd_array_t<std::size_t> gap_idx_{local_box().shape()};
    nd_array_t<std::size_t> add_into_{local_box().shape()};
    nd_array_t<std::size_t> left_{local_box().shape()};
    nd_array_t<std::size_t> cap_{local_box().shape()};
    nd_array_t<std::size_t> cell_size_{local_box().shape()};

    std::size_t total_size = 0;
};


// Span-side counterpart of ParticleTranslationTracker: the same box/local-cell plumbing plus
// the subset of bookkeeping views identical across PerCellSpan / TileSetSpan / PCTileSetSpan
// (gaps_/gap_idx_/add_into_/cap_/left_, size_, box_/ghost_box_/local_ghost_box_).
//
// cell_size_ is deliberately excluded: PerCellSpan has no such member (it tracks per-cell size
// via each cell's own vector), so it isn't identical across all three - it stays declared
// locally on TileSetSpan/PCTileSetSpan instead. Same reasoning as PerCellVector's off_sets_ on
// the VECTOR side above.
template<std::size_t dim>
class ParticleTranslationTrackerSpan
{
public:
    using lobox_t = Box<std::uint32_t, dim>;

    auto size() const { return size_; }

    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }

    auto local_cell(std::array<int, dim> const& icell) const
    {
        return as_local_cell(ghost_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const { return local_cell(icell.toArray()); }

    auto& local_box() const { return local_ghost_box_; }
    auto local_box(Box<int, dim> const& from) const
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

protected:
    // named fields (rather than positional ctor args) so callers can't silently swap two
    // same-typed values (e.g. box vs ghost_box, or the several std::size_t views)
    struct Init
    {
        NdArrayView<dim, Span<std::size_t>> gaps;
        NdArrayView<dim, std::size_t> gap_idx, add_into, cap, left;
        std::size_t size;
        Box<int, dim> box, ghost_box;
        lobox_t local_ghost_box;
    };

    ParticleTranslationTrackerSpan(Init const& init)
        : gaps_{init.gaps}
        , gap_idx_{init.gap_idx}
        , add_into_{init.add_into}
        , cap_{init.cap}
        , left_{init.left}
        , size_{init.size}
        , box_{init.box}
        , ghost_box_{init.ghost_box}
        , local_ghost_box_{init.local_ghost_box}
    {
    }

    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, std::size_t> gap_idx_, add_into_, cap_, left_;
    std::size_t size_;

    Box<int, dim> box_, ghost_box_;
    lobox_t local_ghost_box_;
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_TRANSLATION_TRACKING_HPP */
