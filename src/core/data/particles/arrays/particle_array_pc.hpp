#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP

#include "core/operators.hpp"
#include "core/utilities/span.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"

#if PHARE_HAVE_THRUST
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/arrays/particle_array_soavx.hpp"
#endif



namespace PHARE::core
{

template<typename Particles>
class PerCellSpan
{
protected:
    using SIZE_T = default_span_size_t;

public:
    auto static constexpr alloc_mode = Particles::alloc_mode;
    auto static constexpr dim        = Particles::dimension;
    using This                       = PerCellSpan<Particles>;
    using lobox_t                    = Box<std::uint32_t, dim>;
    using per_cell_particles         = Particles;

private:
    using locell_t = std::array<std::uint32_t, dim>;

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

public:
    auto static constexpr dimension    = dim;
    auto static constexpr storage_mode = StorageMode::SPAN;
    using Particle_t                   = typename ParticleDefaults<dim>::Particle_t;
    using view_t                       = PerCellSpan<Particles>;

    template<typename PerCellArray>
    PerCellSpan(PerCellArray& arr)
        : particles_{resolve(arr)}
        , gaps_{resolve_gaps(arr)}
        , off_sets_{arr.off_sets_}
        , gap_idx_{arr.gap_idx_}
        , add_into_{arr.add_into_}
        , cap_{arr.cap_}
        , left_{arr.left_}
        , p2c_{arr.p2c_.data(), arr.p2c_.size()}
        , size_{arr.size()}
        , box_{arr.box_}
        , ghost_box_{arr.ghost_box()}
        , local_ghost_box_{arr.local_box()}
    {
    }

    auto size() const _PHARE_ALL_FN_ { return size_; }
    auto size(std::size_t const& idx) const _PHARE_ALL_FN_ { return particles_.data()[idx].size(); }

    auto& box() const _PHARE_ALL_FN_ { return box_; }
    auto& ghost_box() const _PHARE_ALL_FN_ { return ghost_box_; }
    auto local_cell(std::array<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return as_local_cell(ghost_box_, icell);
    }
    auto local_cell(Point<int, dim> const& icell) const _PHARE_ALL_FN_
    {
        return local_cell(icell.toArray());
    }

    auto& local_box() const _PHARE_ALL_FN_ { return local_ghost_box_; }

    auto local_box(Box<int, dim> const& from) const _PHARE_ALL_FN_
    {
        return lobox_t{local_cell(from.lower), local_cell(from.upper)};
    }

    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) _PHARE_ALL_FN_ { return particles_(cell); }
    auto& operator()(locell_t const& cell) const _PHARE_ALL_FN_ { return particles_(cell); }

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync() _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain, typename... Args>
    void sync(Args&&... args) _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_add_new() _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left() _PHARE_ALL_FN_;

    // per-cell equivalents of the above, callable from host or device
    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_add_new(locell_t const& bix) _PHARE_ALL_FN_;

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left(locell_t const& bix) _PHARE_ALL_FN_;

    // as above, but also removes an externally registered, ascending-sorted list of
    // leavers in the same descending swap-pop — two separate passes would move particles
    // the other list still indexes
    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync_rm_left(locell_t const& bix, std::size_t const* x_gaps, SIZE_T& x_size,
                      SIZE_T& x_left) _PHARE_ALL_FN_;

    // append src[idx] into cell bix if capacity allows — false when the cell is full
    bool append_from(locell_t const& bix, auto const& src, std::size_t const idx) _PHARE_ALL_FN_
    {
        bool constexpr static GPU = alloc_mode == AllocatorMode::GPU_UNIFIED;

        auto const& cap = cap_(bix);
        auto& nparts    = particles_(bix);
        using Op        = Operators<std::decay_t<decltype(*nparts.size_address())>, true, GPU>;
        auto npidx      = nparts.size();
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

    void static sort(auto from, auto to) _PHARE_ALL_FN_
    {
        [[maybe_unused]] int sorted = 0;
        PHARE_WITH_THRUST({ //
            ++sorted;
            thrust::sort(thrust::seq, from, to);
        })
        //
        PHARE_WITH_THRUST_ELSE({ //
            ++sorted;
            std::sort(from, to);
        }) //

        assert(sorted == 1);
    }

    void clear()
    {
        for (auto const& bix : local_box())
            particles_(bix).clear();
        size_ = 0;
    }

    void reset_p2c(std::array<std::uint32_t, dim>* cells, std::size_t size)
    {
        p2c_ = Span<std::array<std::uint32_t, dim>>{cells, size};
    }


protected:
    NdArrayView<dim, Particles> particles_;
    NdArrayView<dim, Span<std::size_t>> gaps_;
    NdArrayView<dim, SIZE_T> off_sets_, gap_idx_, add_into_, cap_, left_;
    Span<std::array<std::uint32_t, dim>> p2c_;
    std::size_t size_;

    Box<int, dim> box_, ghost_box_;
    lobox_t local_ghost_box_;


}; // PerCellSpan


template<typename Particles>
class PerCellVector
{
    using This = PerCellVector<Particles>;

protected:
    using SIZE_T = default_span_size_t;

private:
    template<typename>
    friend class PerCellSpan;

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

    template<typename T>
    using vec_helper = PHARE::Vector<T, alloc_mode, 1>;

    using size_t_vec_helper   = vec_helper<std::size_t>;
    using size_t_vector       = typename size_t_vec_helper::vector_t;
    using particle_vec_helper = vec_helper<Particle_t>;
    using particle_vector     = typename particle_vec_helper::vector_t;

    PerCellVector(box_t const& box = {}, std::size_t ghost_cells = 0)
        : ghost_cells_{ghost_cells}
        , box_{box}
        , ghost_box_{grow(box, ghost_cells)}
    {
        cell_size_.zero();
        gap_idx_.zero();
        add_into_.zero();
        left_.zero();
        cap_.zero();

        reset_views();
    }

    PerCellVector(PerCellVector const& from)            = default;
    PerCellVector(PerCellVector&& from)                 = default;
    PerCellVector& operator=(PerCellVector&& from)      = default;
    PerCellVector& operator=(PerCellVector const& from) = default;

    template<auto type = ParticleType::Domain>
    auto& reserve_ppc(std::size_t const& ppc);

    void reserve(std::size_t const&) {} // per-cell storage pre-allocated; hint is a no-op

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
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            PHARE_WITH_MKN_GPU(return mkn::gpu::as_super(v));
        }
        else
        {
            return v;
        }
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


    auto& box() const { return box_; }
    auto& ghost_box() const { return ghost_box_; }


    auto& operator()() const { return particles_; }
    auto& operator()(locell_t const& cell) { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) { return particles_.data() + cell; }
    auto& operator()(locell_t const& cell) const { return particles_(cell); }
    auto& operator()(std::uint32_t const& cell) const { return particles_.data() + cell; }

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

    template<std::uint8_t PHASE = 0, auto type = ParticleType::Domain>
    void sync(); // after move

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

    template<auto type>
    void trim();



    void clear()
    {
        for (auto const& bix : local_box())
            particles_(bix).clear();
        total_size = 0;
    }


    auto& insert(PerCellVector const& src);


    void static resize(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (any_in(layout_mode, LayoutMode::SoA))
            apply(ps.as_tuple(), [&](auto& v) { resize(v, s, copy); });
        else
            resize(ps.particles_, s, copy);
    }

    template<typename V>
    void static resize(V& v, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(mkn::gpu::resize(v, s, copy));

        else
            v.resize(s);
    }

    void static reserve(Particles& ps, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (any_in(layout_mode, LayoutMode::SoA))
            apply(ps.as_tuple(), [&](auto& v) { reserve(v, s, copy); });
        else
            reserve(ps.particles_, s, copy);
    }

    template<typename V>
    void static reserve(V& v, std::size_t const& s, bool const& copy = true)
    {
        if constexpr (CompileOptions::WithMknGpu and alloc_mode == AllocatorMode::GPU_UNIFIED)
            PHARE_WITH_MKN_GPU(mkn::gpu::reserve(v, s, copy));

        else
            v.reserve(s);
    }

    template<typename View_t>
    void reset_p2c(View_t& view)
    {
        view.reset_p2c(p2c_.data(), p2c_.size());
    }

protected:
    bool static constexpr c_order = true;

    std::size_t ghost_cells_;

    Box<int, dim> box_, ghost_box_;

    NdArrayVector<dim, Particles, c_order, alloc_mode> particles_{local_box().shape()};

    // only used for GPU

    NdArrayVector<dim, PSpan_t, c_order, alloc_mode> particles_views_{local_box().shape()};

    typename vec_helper<std::array<std::uint32_t, dim>>::vector_t p2c_;
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> off_sets_{local_box().shape()};

    NdArrayVector<dim, size_t_vector, c_order, alloc_mode> gaps_{local_box().shape()};
    NdArrayVector<dim, Span<std::size_t>, c_order, alloc_mode> gap_views_{local_box().shape()};

    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> gap_idx_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> add_into_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> left_{local_box().shape()};
    NdArrayVector<dim, SIZE_T, c_order, alloc_mode> cap_{local_box().shape()};
    NdArrayVector<dim, std::size_t, c_order, alloc_mode> cell_size_{local_box().shape()};

    std::size_t total_size = 0;


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
        sync<2>();

    return *this;
}

template<typename Particles>
template<auto type>
auto& PerCellVector<Particles>::reserve_ppc(std::size_t const& ppc)
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    std::size_t const additional = ppc < 50 ? 10 : ppc * .2; // 20% overallocate
    std::size_t const buffered   = ppc + additional;

    if constexpr (type == ParticleType::Domain)
    {
        on_domain([&](auto const& bix) { particles_(bix).reserve(buffered); });
        on_ghost_box([&](auto const& bix) { reserve(gaps_(bix), additional); });
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(additional); });
    }

    if constexpr (type == ParticleType::Ghost)
    {
        on_ghost_layer_plus_2_domain([&](auto const& bix) {
            particles_(bix).reserve(additional);
            reserve(gaps_(bix), additional);
        });
        on_ghost_layer([&](auto const& bix) { particles_(bix).reserve(buffered); });
    }

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
    cap_(bix)      = real.capacity();
    add_into_(bix) = 0;
}


template<typename Particles>
template<auto type>
void PerCellVector<Particles>::trim()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);

    if constexpr (type == ParticleType::Ghost)
        for (auto const& ghost_layer_box : local_box().remove(local_box(box_)))
            for (auto const& bix : ghost_layer_box)
            {
                cell_size_(bix) = 0;
                particles_(bix).clear();
                gaps_(bix).clear();
            }

    if constexpr (type == ParticleType::Domain)
    {
        throw std::runtime_error("NO");
    }
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
template<std::uint8_t PHASE, auto type>
void PerCellVector<Particles>::sync()
{
    static_assert(std::is_same_v<decltype(type), ParticleType>);
    static_assert(type != ParticleType::All);

    PHARE_LOG_SCOPE(1, "PerCellVector::sync");
    // PHARE_LOG_LINE_STR("sync " << static_cast<std::uint32_t>(PHASE) << " "
    //                            << magic_enum::enum_name(type));

    auto const lbox = local_box();

    total_size = 0;
    for (auto const& bix : lbox)
        total_size += (cell_size_(bix) = particles_(bix).size());

    PHARE_LOG_SCOPE(1, "PerCellVector::sync::reset");

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

    {
        PHARE_LOG_SCOPE(1, "PerCellVector::sync::reset_views");
        reset_views();
    }
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

    template<typename T>
    struct iterator_impl;

    template<typename T, typename... Args>
    auto static it(T* t, Args&&... args)
    {
        if constexpr (storage_mode == StorageMode::SPAN)
            return iterator_impl<T>{*t, args...};
        else
            return iterator_impl<T&>{*t, args...};
    }
    auto begin() const _PHARE_ALL_FN_ { return it(this); }
    auto begin() _PHARE_ALL_FN_ { return it(this); }
    auto end() const _PHARE_ALL_FN_ { return it(this, true); }
    auto end() _PHARE_ALL_FN_ { return it(this, true); }

    auto data() const _PHARE_ALL_FN_ { return particles_.data(); }
    auto data() _PHARE_ALL_FN_ { return particles_.data(); }



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
    auto& move_check(auto const& pt, std::size_t const& idx, auto& particle) _PHARE_ALL_FN_
    {
        static_assert(any_in(particle_type, ParticleType::Domain, ParticleType::LevelGhost));

        auto const& newcell = particle.iCell();
        if (array_equals(newcell, pt.icell))
            return *this; // old cell == new cell, no change required

        bool constexpr static ATOMIC = true;
        bool constexpr static GPU    = alloc_mode == AllocatorMode::GPU_UNIFIED;

        using Op = Operators<typename Super::SIZE_T, ATOMIC, GPU>;

        auto const old_lcl_cell = Super::local_cell(pt.icell);
        Super::gaps_(old_lcl_cell)[Op{Super::gap_idx_(old_lcl_cell)}.increment_return_old()] = idx;

        if (isIn(newcell, Super::ghost_box()))
            Op{Super::add_into_(Super::local_cell(newcell))}.increment_return_old();

        return *this;
    }

    void print() const {}
    void check() const {}

    auto max_size() const
    {
        return max_from(this->particles_,
                        [](auto const& v, auto const& i) { return v.data()[i].size(); });
    }


}; // PerCellParticles<Super>




template<typename OuterSuper>
template<typename T>
struct PerCellParticles<OuterSuper>::iterator_impl
{
    auto static constexpr dimension = OuterSuper::dimension;

    using outer_type        = std::decay_t<T>;
    using difference_type   = std::size_t;
    using iterator_category = std::forward_iterator_tag;
    using Particle_t        = typename OuterSuper::Particle_t;
    using value_type        = Particle_t;
    using pointer           = Particle_t*;
    using reference         = Particle_t&;

    iterator_impl(T& particles_, bool end = false)
        : particles{particles_}
    {
        if (end)
        {
            l = particles.particles_.size() - 1;
            i = particles.size(l);
        }
        else
            inc(); // find first particle
    }
    iterator_impl(iterator_impl&& that)                 = default;
    iterator_impl(iterator_impl const& that)            = default;
    iterator_impl& operator=(iterator_impl&& that)      = default;
    iterator_impl& operator=(iterator_impl const& that) = default;

    auto end() { return iterator_impl{particles, true}; }
    auto& inc()
    {
        i = 0;
        for (; l < particles.particles_.size(); ++l)
            if (*this == this->end() or particles.size(l) > 0)
                break;
        return *this;
    }

    iterator_impl& operator++()
    {
        auto last = end();
        if (*this == last)
            return *this;
        ++i;
        if (*this == last)
            return *this;
        if (i == particles.size(l))
        {
            ++l;
            return inc();
        }
        return *this;
    }
    auto operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    auto operator==(iterator_impl const& that) const { return l == that.l and i == that.i; }
    auto operator!=(iterator_impl const& that) const { return !(*this == that); }
    auto operator<(iterator_impl const& that) const
    {
        if (l < that.l)
            return true;
        if (l == that.l)
            return i < that.i;
        return false;
    }

    auto& charge() _PHARE_ALL_FN_ { return particles.data()[l][i].charge(); }
    auto& delta() _PHARE_ALL_FN_ { return particles.data()[l][i].delta(); }
    auto& iCell() _PHARE_ALL_FN_ { return particles.data()[l][i].iCell(); }
    auto& weight() _PHARE_ALL_FN_ { return particles.data()[l][i].weight(); }
    auto& v() _PHARE_ALL_FN_ { return particles.data()[l][i].v(); }

    auto& charge() const _PHARE_ALL_FN_ { return particles.data()[l][i].charge(); }
    auto& delta() const _PHARE_ALL_FN_ { return particles.data()[l][i].delta(); }
    auto& iCell() const _PHARE_ALL_FN_ { return particles.data()[l][i].iCell(); }
    auto& weight() const _PHARE_ALL_FN_ { return particles.data()[l][i].weight(); }
    auto& v() const _PHARE_ALL_FN_ { return particles.data()[l][i].v(); }

    auto& operator*() _PHARE_ALL_FN_ { return particles.data()[l][i]; }
    auto& operator*() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    auto copy() const _PHARE_ALL_FN_ { return particles.data()[l][i]; }

    T particles;
    std::size_t l = 0, i = 0;
};




template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync() _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(1, "PerCellSpan::sync()");

    auto view = *this;
    PHARE_WITH_MKN_GPU({
        auto const& gabox = local_ghost_box_;
        mkn::gpu::GDLauncher{gabox.size()}(
            [=] _PHARE_ALL_FN_() mutable { view.template sync_add_new<PHASE, type>(); });
        mkn::gpu::GDLauncher{gabox.size()}(
            [=] _PHARE_ALL_FN_() mutable { view.template sync_rm_left<PHASE, type>(); });
    })
}

template<typename Particles>
template<std::uint8_t PHASE, auto type, typename... Args>
void PerCellSpan<Particles>::sync(Args&&... args) _PHARE_ALL_FN_
{
    PHARE_LOG_SCOPE(1, "PerCellSpan::sync(stream)");

    auto view = *this;

    if constexpr (alloc_mode == AllocatorMode::CPU)
    {
        for (auto const& bix : local_ghost_box_)
            sync_add_new<PHASE, type>(bix.toArray());

        for (auto const& bix : local_ghost_box_)
            sync_rm_left<PHASE, type>(bix.toArray());
    }
    else if (alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        PHARE_WITH_MKN_GPU({
            auto const& [stream] = std::forward_as_tuple(args...);
            auto const& gabox    = local_ghost_box_;
            mkn::gpu::GDLauncher<true>{gabox.size()}.stream(stream, [=] _PHARE_ALL_FN_() mutable {
                view.template sync_add_new<PHASE, type>();
            });
            mkn::gpu::GDLauncher<true>{gabox.size()}.stream(stream, [=] _PHARE_ALL_FN_() mutable {
                view.template sync_rm_left<PHASE, type>();
            });
        })
    }
}


template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_add_new() _PHARE_ALL_FN_
{
#if PHARE_HAVE_MKN_GPU
    sync_add_new<PHASE, type>((*(local_ghost_box_.begin() + mkn::gpu::idx())).toArray());
#endif // PHARE_HAVE_MKN_GPU
}

template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_add_new(locell_t const& bix) _PHARE_ALL_FN_
{
    auto const& n_gaps = gap_idx_(bix);
    {
        auto& gaps = gaps_(bix);
        sort(gaps.data(), gaps.data() + n_gaps /*, std::greater<>()*/);
    }
    auto& real       = particles_(bix);
    auto& left       = left_(bix);
    auto const& gaps = gaps_(bix);
    for (std::size_t i = 0; i < n_gaps; ++i)
    {
        auto const& gidx = gaps[n_gaps - (1 + i)];
        if (!append_from(local_cell(real.iCell(gidx)), real, gidx))
            return;
        ++left;
    }
}


template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_rm_left() _PHARE_ALL_FN_
{
#if PHARE_HAVE_MKN_GPU
    sync_rm_left<PHASE, type>((*(local_ghost_box_.begin() + mkn::gpu::idx())).toArray());
#endif // PHARE_HAVE_MKN_GPU
}

template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_rm_left(locell_t const& bix) _PHARE_ALL_FN_
{
    add_into_(bix) -= left_(bix);
    SIZE_T x_size = 0, x_left = 0;
    sync_rm_left<PHASE, type>(bix, nullptr, x_size, x_left);
}

template<typename Particles>
template<std::uint8_t PHASE, auto type>
void PerCellSpan<Particles>::sync_rm_left(locell_t const& bix, std::size_t const* x_gaps,
                                          SIZE_T& x_size, SIZE_T& x_left) _PHARE_ALL_FN_
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


template<auto type>
void sync_moved_pc(auto& particles, auto&&... args)
{
    particles.template sync_moved<type>();   // realloc
    particles.template sync_add_new<type>(); // movers appended to their new cells
    particles.template sync_rm_left<type>(); // leavers compacted out of their old cells

    particles.template sync<0, type>(args...); // finalize
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_PerCell_HPP */
