#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP

#include "core/operators.hpp"
#include "core/utilities/span.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_translation_tracking.hpp"

namespace PHARE::core
{

template<typename Particles>
class PerCellSpan : public ParticleTranslationTrackerSpan<Particles::dimension>
{
    using Base = ParticleTranslationTrackerSpan<Particles::dimension>;

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
    using This                       = PerCellSpan<Particles>;
    using lobox_t                    = Box<std::uint32_t, dim>;
    using per_cell_particles         = Particles;
    using locell_t                   = std::array<std::uint32_t, dim>;
    using Base::box;
    using Base::ghost_box;
    using Base::local_box;
    using Base::local_cell;
    using Base::size;

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;
    using view_t                       = PerCellSpan<Particles>;

    template<typename PerCellArray>
    PerCellSpan(PerCellArray& arr);

    auto size(std::size_t const& idx) const { return particles_.data()[idx].size(); }

    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return particles_(cell); }
    auto& operator()(locell_t const& cell) const { return particles_(cell); }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain, typename... Args>
    void sync(Args&&... args);

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_add_new(locell_t const& bix);

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left(locell_t const& bix);

    // as above, but also removes an externally registered, ascending-sorted list of
    // leavers in the same descending swap-pop — two separate passes would move particles
    // the other list still indexes
    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left(locell_t const& bix, std::size_t const* x_gaps, std::size_t& x_size,
                      std::size_t& x_left);

    // append src[idx] into cell bix if capacity allows — false when the cell is full
    bool append_from(locell_t const& bix, auto const& src, std::size_t const idx);

    void static sort(auto from, auto to) { std::sort(from, to); }

    void clear()
    {
        for (auto const& bix : local_box())
            particles_(bix).clear();
        size_ = 0;
    }

protected:
    template<typename PerCellArray>
    auto resolve(PerCellArray& arr)
    {
        if constexpr (PerCellArray::storage_mode == StorageMode::SPAN)
            return arr.particles_;
        else
        {
            arr.check();
            return *arr.particles_views_;
        }
    }

    template<typename PerCellArray>
    auto resolve_gaps(PerCellArray& arr)
    {
        if constexpr (PerCellArray::storage_mode == StorageMode::SPAN)
            return arr.gaps_;
        else
        {
            arr.check();
            return *arr.gap_views_;
        }
    }

    NdArrayView<dim, Particles> particles_;
    NdArrayView<dim, std::size_t> off_sets_;

}; // PerCellSpan


template<typename Particles>
class PerCellVector : public ParticleTranslationTracker<Particles::dimension, Particles::alloc_mode>
{
    template<typename>
    friend class PerCellSpan;

    bool static constexpr c_order = true;
    using This                    = PerCellVector<Particles>;
    using Base = ParticleTranslationTracker<Particles::dimension, Particles::alloc_mode>;

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

public:
    auto static constexpr dim          = Particles::dimension;
    auto static constexpr alloc_mode   = Particles::alloc_mode;
    auto static constexpr layout_mode  = Particles::layout_mode;
    auto static constexpr storage_mode = StorageMode::VECTOR;
    auto static constexpr dimension    = dim;
    using box_t                        = Box<int, dim>;
    using lobox_t                      = Box<std::uint32_t, dim>;
    using locell_t                     = std::array<std::uint32_t, dim>;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;
    using value_type                   = Particle_t;
    using PSpan_t                      = typename Particles::view_t;
    using view_t                       = PerCellSpan<PSpan_t>;
    using per_cell_particles           = Particles;
    using per_tile_particles           = Particles; // MultiBoris compatibility

    using size_t_vector   = std::vector<std::size_t>;
    using particle_vector = std::vector<Particle_t>;

    using Base::box;
    using Base::ghost_box;
    using Base::local_box;
    using Base::local_cell;
    using Base::total_size;
    using Base::zero_bookkeeping;

    PerCellVector(box_t const& box = {}, std::size_t ghost_cells = 0)
        : Base{box, ghost_cells}
    {
        reset_views();
    }

    PerCellVector(PerCellVector const& from)            = default;
    PerCellVector(PerCellVector&& from)                 = default;
    PerCellVector& operator=(PerCellVector&& from)      = default;
    PerCellVector& operator=(PerCellVector const& from) = default;

    auto size() const { return total_size; }
    auto size(std::array<std::uint32_t, dim> const& icell) const { return cell_size_(icell); }
    auto size(std::size_t const& idx) const { return cell_size_.data()[idx]; }

    void _inc(locell_t const& locell)
    {
        ++total_size;
        ++cell_size_(locell);
    }

    template<bool inc_ = true>
    void emplace_back(Particle_t const& p)
    {
        auto const locell = local_cell(p.iCell());
        particles_(locell).emplace_back(p);
        if constexpr (inc_)
            _inc(locell);
    }

    template<bool inc_ = true>
    void emplace_back(Particles& dst, Particles const& src, std::size_t const& idx)
    {
        dst.emplace_back(src, idx);
    }


    template<typename... Args>
    void emplace_back(double const weight, Args&&... args)
    {
        this->emplace_back(Particle_t{weight, args...});
    }

    template<typename V>
    static auto& get_vec(V& v)
    {
        return v;
    }

    void push_back(Particle_t&& p) { emplace_back(p); }
    void push_back(Particle_t const& p) { emplace_back(p); }

    void reset_views()
    {
        particles_views_ = generate_from<alloc_mode>(
            [&](auto const i) {
                return PSpan_t{*(particles_.data() + i), *(cell_size_.data() + i)};
            },
            particles_);
        gap_views_ = generate_from<alloc_mode>(
            [&](auto const i) { return make_span(*(gaps_.data() + i)); }, gaps_);
    }


    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

    // refreshes recount/gap-sizing/views after particles were added without going through
    // move_check (a raw append) - see on_moved() for the post-move_check pipeline.
    template<auto type = ParticleType::Domain>
    void on_appended();

    // applies move_check-registered moves: realloc/resize, append movers, compact leavers,
    // then on_appended() to finalize.
    template<auto type>
    void on_moved();

    template<auto type>
    void sync_moved(); // realloc for particles registered as incoming via add_into_

    // realloc + resize one cell ahead of span-side copies; in/out are externally
    // registered arrivals/departures on top of this cell's own counters; records the
    // pre-move size so freshly reset views start the copy step from it
    template<auto type>
    void sync_moved(locell_t const& bix, std::size_t const in, std::size_t const out);

    template<auto type>
    void sync_add_new(); // append registered movers into their new cells

    template<auto type>
    void sync_rm_left(); // swap-delete registered leavers, reset gap/add counters

    void clear();

    auto& insert(PerCellVector const& src);

    void static resize(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        resize(ps.particles_, s, copy);
    }

    template<typename V>
    void static resize(V& v, std::size_t const& s, bool const& copy = true)
    {
        v.resize(s);
    }

    void static reserve(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        reserve(ps.particles_, s, copy);
    }

    template<typename V>
    void static reserve(V& v, std::size_t const& s, bool const& copy = true)
    {
        v.reserve(s);
    }

protected:
    static void on_box(auto&& box, auto&& fn)
    {
        for (auto const& bix : box)
            fn(bix);
    };

    static void on_box_list(auto&& boxlist, auto&& fn)
    {
        for (auto const& box : boxlist)
            on_box(box, fn);
    };

    void on_domain(auto&& fn) const { on_box(local_box(box()), fn); };
    void on_ghost_box(auto&& fn) { on_box(local_box(), fn); };
    void on_ghost_layer(auto&& fn) const { on_box_list(local_box().remove(local_box(box())), fn); };
    void on_ghost_layer_plus_2_domain(auto&& fn) const
    {
        on_box_list(local_box().remove(shrink(local_box(box()), 2)), fn);
    };


    NdArrayVector<dim, Particles, c_order, alloc_mode> particles_{local_box().shape()};
    NdArrayVector<dim, PSpan_t, c_order, alloc_mode> particles_views_{local_box().shape()};
    NdArrayVector<dim, std::size_t, c_order, alloc_mode> off_sets_{local_box().shape()};

}; // PerCellVector<Particles>


template<typename Particles>
auto& PerCellVector<Particles>::insert(PerCellVector const& src)
{
    std::size_t added = 0;
    for (auto const& bix : local_box(box()))
    {
        auto& from = src(bix);
        auto& to   = (*this)(bix);
        added += from.size();
        to.reserve(to.size() + from.size());
        std::copy(from.begin(), from.end(), std::back_inserter(to));
    }
    if (added)
        on_appended();

    return *this;
}


template<typename Particles>
template<auto type>
void PerCellVector<Particles>::sync_moved()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    for (auto const& bix : local_box())
    {
        auto& real      = particles_(bix);
        cell_size_(bix) = real.size();
        resize(real, real.size() + add_into_(bix)); // arrivals assigned in sync_add_new
    }
}

template<typename Particles>
template<auto type>
void PerCellVector<Particles>::sync_moved(locell_t const& bix, std::size_t const in,
                                          std::size_t const out)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    auto& real          = particles_(bix);
    auto const old_size = real.size();
    auto const nu       = add_into_(bix) + in;
    cell_size_(bix)     = old_size;
    reserve(real, old_size + nu); // copies transiently exceed the final size before rm
    resize(real, old_size + nu - (gap_idx_(bix) + out));
    cap_(bix) = real.capacity();
    // add_into_(bix) is deliberately left as-is here (not zeroed): sync_add_new still
    // needs it as the precise "is an arrival actually expected at this cell" signal
    // (move_check only increments it when it intends the particle to land here — a
    // departure register_particle_ejection registered never does) — a plain isIn(ghost_box())
    // re-check can't tell "stays local" from "crossed into a sibling tile's own
    // domain that happens to overlap this tile's ghost reach" apart, since both are
    // geometrically within ghost_box(). sync_rm_left resets it to 0 once consumed.
}


// add incoming particles - move_check only registered the move, the particle's own
// iCell() was already updated by the caller before registering, so it tells us where
// each registered departure from this cell should actually land
template<typename Particles>
template<auto type>
void PerCellVector<Particles>::sync_add_new()
{
    // arrivals are assigned into the slots opened by sync_moved's resize, cursored by
    // left_ — never emplaced, so per-cell sizes stay stable for the rm pass
    for (auto const& bix : local_box())
    {
        auto& real            = particles_(bix);
        auto const& gaps      = gaps_(bix);
        auto const& gaps_size = gap_idx_(bix);
        for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
        {
            auto const& idx   = gaps[gidx];
            auto const& icell = real.iCell(idx);
            if (not isIn(icell, ghost_box()))
                continue; // ghost-box leaver (level ghost): no destination bucket,
                          // the rm pass deletes it from its old cell
            auto const newcell = local_cell(icell);
            particles_(newcell).assign(real, idx, cell_size_(newcell) + left_(newcell)++);
        }
    }
}

// delete outgoing particles from their old cell (descending order keeps indices valid)
template<typename Particles>
template<auto type>
void PerCellVector<Particles>::sync_rm_left()
{
    for (auto const& bix : local_box())
    {
        auto const& gaps_size = gap_idx_(bix);
        {
            auto& gaps = gaps_(bix);
            std::sort(gaps.begin(), gaps.begin() + gaps_size, std::greater<>());
        }
        auto& real       = particles_(bix);
        auto const& gaps = gaps_(bix);
        for (std::size_t gidx = 0; gidx < gaps_size; ++gidx)
        {
            auto const& idx = gaps[gidx];
            real.assign(real.size() - 1, idx);
            real.pop_back();
        }
        gap_idx_(bix)  = 0;
        add_into_(bix) = 0;
        left_(bix)     = 0;
    }
}

template<typename Particles>
template<auto type>
void PerCellVector<Particles>::on_appended()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

    PHARE_LOG_SCOPE(3, "PerCellVector::on_appended");

    total_size = 0;
    for (auto const& bix : local_box())
        total_size += (cell_size_(bix) = particles_(bix).size());

    auto const per_cell = [&](auto& bix) {
        auto const& cs  = cell_size_(bix);
        auto const& cap = particles_(bix).capacity();
        auto& gaps      = gaps_(bix);
        if (gaps.size() < cs)
        {
            reserve(gaps, cap, false);
            resize(gaps, cs, false);
        }
        cap_(bix) = cap;
    };

    if constexpr (type == ParticleType::Domain)
        on_ghost_box(per_cell);
    else
        on_ghost_layer_plus_2_domain(per_cell);

    reset_views();
}

template<typename Particles>
void PerCellVector<Particles>::clear()
{
    for (auto const& bix : local_box())
        particles_(bix).clear();
    total_size = 0;

    // iteration/size() go through cell_size_, not each cell vector's own size(), so
    // leaving it stale here makes clear() a no-op from any caller's perspective -
    // reset every bit of per-cell bookkeeping the same way the constructor does.
    zero_bookkeeping();
    reset_views();
}


template<typename Super_>
struct PerCellParticles : public Super_
{
    using Super              = Super_;
    using This               = PerCellParticles<Super>;
    using Particle_t         = typename Super::Particle_t;
    using per_cell_particles = typename Super::per_cell_particles;
    using view_t             = PerCellParticles<typename Super::view_t>;

    auto static constexpr alloc_mode   = Super::alloc_mode;
    auto static constexpr dimension    = Super::dimension;
    auto static constexpr storage_mode = Super::storage_mode;
    auto static constexpr size_of_particle() { return sizeof(Particle_t); }

    using Super::local_box;
    using Super::particles_;
    using Super::size;

    template<typename... Args>
    PerCellParticles(Args&&... args)
        : Super{std::forward<Args>(args)...}
    {
    }

    PerCellParticles(PerCellParticles const& from)            = default;
    PerCellParticles(PerCellParticles&& from)                 = default;
    PerCellParticles& operator=(PerCellParticles&& from)      = default;
    PerCellParticles& operator=(PerCellParticles const& from) = default;

    auto data() const { return particles_.data(); }
    auto data() { return particles_.data(); }

    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](Box<std::uint32_t, dimension> const& local) const
    {
        std::vector<Particle_t> out;
        out.reserve(sum_from(local, [&](auto const& b) { return particles_(b.toArray()).size(); }));
        for (auto const& b : local)
            std::copy(particles_(b).begin(), particles_(b).end(), std::back_inserter(out));
        return out;
    }
    template<auto S = storage_mode, typename = std::enable_if_t<S == StorageMode::VECTOR>>
    auto operator[](Box<int, dimension> const& amr) const
    {
        return (*this)[Super::local_box(amr)];
    }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }

    // the particle already carries its NEW cell (set by the caller); pt carries the OLD
    // cell so the departure can be registered against it — mirrors the TS/PCTS move_check
    template<auto particle_type>
    auto& move_check(auto const& pt, std::size_t const& idx, auto& particle);

    // the new cell isn't ours, even if still within our ghost_box() reach (a
    // neighbour's domain can overlap it) — register the departure only, no arrival
    void register_particle_ejection(auto const& pt, std::size_t const& idx)
    {
        bool constexpr static ATOMIC = true;
        using Op                     = Operators<std::size_t, ATOMIC>;

        auto const old_lcl_cell = Super::local_cell(pt.icell);
        Super::gaps_(old_lcl_cell)[Op{Super::gap_idx_(old_lcl_cell)}.increment_return_old()] = idx;
    }

    void print() const {}
    void check() const {}

    auto max_size() const
    {
        return max_from(this->particles_,
                        [](auto const& v, auto const& i) { return v.data()[i].size(); });
    }
}; // PerCellParticles<Super>


// the particle already carries its NEW cell (set by the caller); pt carries the OLD
// cell so the departure can be registered against it — mirrors the TS/PCTS move_check
template<typename Super_>
template<auto particle_type>
auto& PerCellParticles<Super_>::move_check(auto const& pt, std::size_t const& idx, auto& particle)
{
    static_assert(any_in(particle_type, ParticleType::Domain, ParticleType::LevelGhost));
    bool constexpr static ATOMIC = true;
    using Op                     = Operators<std::size_t, ATOMIC>;

    auto const& newcell = particle.iCell();
    if (array_equals(newcell, pt.icell))
        return *this; // old cell == new cell, no change required

    auto const old_lcl_cell = Super::local_cell(pt.icell);
    Super::gaps_(old_lcl_cell)[Op{Super::gap_idx_(old_lcl_cell)}.increment_return_old()] = idx;
    if (isIn(newcell, Super::ghost_box()))
        Op{Super::add_into_(Super::local_cell(newcell))}.increment_return_old();

    return *this;
}


template<typename Particles>
template<typename PerCellArray>
PerCellSpan<Particles>::PerCellSpan(PerCellArray& arr)
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
    , off_sets_{arr.off_sets_}
{
}

template<typename Particles>
bool PerCellSpan<Particles>::append_from(locell_t const& bix, auto const& src,
                                         std::size_t const idx)
{
    bool constexpr static atomic = true;
    auto const& cap              = cap_(bix);
    auto& nparts                 = particles_(bix);
    using Op   = Operators<std::decay_t<decltype(*nparts.size_address())>, atomic>;
    auto npidx = nparts.size();
    while (true)
    {
        if (npidx >= cap)
            return false;
        auto const old = Op::compare_and_swap(nparts.size_address(), npidx, npidx + 1);
        if (npidx != old)
        {
            ++npidx;
            continue;
        }
        break;
    }
    nparts.assign(src, idx, npidx);
    return true;
}

template<typename Particles>
template<std::uint8_t PHASE, auto type, typename... Args>
void PerCellSpan<Particles>::sync(Args&&... args)
{
    for (auto const& bix : local_box())
        sync_add_new<PHASE, type>(bix.toArray());

    for (auto const& bix : local_box())
        sync_rm_left<PHASE, type>(bix.toArray());
}


template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_add_new(locell_t const& bix)
{
    auto const& n_gaps = gap_idx_(bix);
    auto& gaps         = gaps_(bix);
    sort(gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
    auto& real = particles_(bix);
    auto& left = left_(bix);
    for (std::size_t i = 0; i < n_gaps; ++i)
    {
        auto const& gidx    = gaps[n_gaps - (1 + i)];
        auto const& newcell = real.iCell(gidx);
        if (not isIn(newcell, ghost_box()))
        { // truly outside this tile's own tracked storage: no destination bucket to
          // append into, but still counts toward sync_rm_left's removal count below
          // (left drives how many of this cell's registered departures actually get
          // swap-popped out)
            ++left;
            continue;
        }
        auto const newcell_local = local_cell(newcell);
        auto& expected           = add_into_(newcell_local);
        if (expected == 0)
        { // still geometrically within this tile's own ghost_box() reach, but no
          // arrival was actually registered for it there — move_check only
          // increments add_into_ at the destination when it intends the particle to
          // land there; register_particle_ejection (a departure crossing into a SIBLING
          // tile's own domain, which this tile's ghost_box() alone can't
          // distinguish from staying local, since a tile's domain overlaps its
          // neighbours' grown ghost box) never does. Trust that signal over blind
          // geometry: drop it here, exactly like a true ghost-box leaver above.
            ++left;
            continue;
        }
        --expected;
        if (!append_from(newcell_local, real, gidx))
            return;
        ++left;
    }
}


template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_rm_left(locell_t const& bix)
{
    // add_into_(bix) is now precisely decremented, once per accepted arrival, right
    // in sync_add_new (see there) — so by the time every tile's own sync_add_new
    // pass has run, it should already be back to 0 for every cell. Reset it
    // unconditionally here (rather than the old `-= left_(bix)`, which conflated
    // bix's arrival count with a same-named but differently-scoped departure count)
    // so the next cycle starts clean regardless.
    add_into_(bix)     = 0;
    std::size_t x_size = 0, x_left = 0;
    sync_rm_left<PHASE, type>(bix, nullptr, x_size, x_left);
}

template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_rm_left(locell_t const& bix, std::size_t const* x_gaps,
                                          std::size_t& x_size, std::size_t& x_left)
{
    auto const& gaps = gaps_(bix);
    auto& real       = particles_(bix);
    auto& left       = left_(bix);
    auto& gaps_size  = gap_idx_(bix);
    while (left or x_left)
    {
        bool const own   = !x_left or (left and gaps[gaps_size - 1] > x_gaps[x_size - 1]);
        auto const& pidx = own ? gaps[gaps_size - 1] : x_gaps[x_size - 1];
        real.assign(real.size() - 1, pidx); // real[pidx] = real[real.size() - 1];
        real.pop_back();
        if (own)
        {
            --gaps_size;
            --left;
        }
        else
        {
            --x_size;
            --x_left;
        }
    }
}

template<typename Particles>
template<auto type>
void PerCellVector<Particles>::on_moved()
{
    sync_moved<type>();   // realloc
    sync_add_new<type>(); // movers appended to their new cells
    sync_rm_left<type>(); // leavers compacted out of their old cells
    on_appended<type>();  // finalize
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP */
