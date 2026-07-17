// IWYU pragma: private, include "core/numerics/pusher/multi_boris.hpp"

#ifndef PHARE_CORE_PUSHER_BORIS_MKN_MULTI_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_MKN_MULTI_BORIS_HPP

#include "core/data/electromag/electromag.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/utilities/span.hpp"


namespace PHARE::core
{

enum class MultiBorisMode : std::uint16_t { REF = 0, COPY };

struct MultiBorisOptions
{
    bool use_main_thread = false; // true for perf
};

// Primary (incomplete) template — specialisations provide the implementation per
// (LayoutMode, AllocatorMode) pair.  Add a new #include below for each new layout.
template<LayoutMode layout, AllocatorMode alloc, typename GridLayout, typename Particles,
         typename Electromag, typename Interpolator>
struct MultiBorisPusherImpl;


// ── MultiBoris state struct for AoSTS / ModelAccessor ─────────────────────────────────

template<typename ModelAccessor, auto _opts = MultiBorisOptions{}>
struct MultiBoris
{
    using Model_t      = ModelAccessor::Model_t;
    using GridLayout_t = ModelAccessor::GridLayout_t;

    static constexpr auto opts = _opts;
    static constexpr auto dim  = GridLayout_t::dimension;


    using ParticleArray_t    = Model_t::particle_array_type;
    using Electromag_t       = Model_t::electromag_type;
    using Vecfield_t         = Electromag_t::vecfield_type;
    using Field_t            = Vecfield_t::field_type;
    using ParticleArray_v    = ParticleArray_t::view_t;
    using Box_t              = Box<int, dim>;
    using Particles_ptrs     = std::vector<ParticleArray_t*>;
    using StreamLauncher     = gpu::ThreadedStreamLauncher<ModelAccessor>;
    using PerTileParticles_t = ParticleArray_t::per_tile_particles;

    MultiBoris(double const dt_, ModelAccessor& _accessor, std::function<void(int)> fn_ = {})
        : dt{dt_}
        , accessor{_accessor}
        , fn{fn_}
    {
    }

    void reset() {}

#if PHARE_HAVE_MKN_GPU
    using GpuBoxAlloc_t = mkn::gpu::ManagedAllocator<Box_t>;
#else
    using GpuBoxAlloc_t = std::allocator<Box_t>;
#endif
    using GpuBoxSpanSet_t = SpanSet<Box_t, default_span_size_t, GpuBoxAlloc_t>;

    double const dt;
    ModelAccessor& accessor;
    std::function<void(int)> fn;
    StreamLauncher streamer{accessor, opts.use_main_thread ? 0 : 1};
    GpuBoxSpanSet_t gpu_nlgb;

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](double const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }
};


// ── Per-tile functors (shared between CPU and GPU specialisations) ─────────────────────

template<auto particle_type, auto boris_mode, typename MultiBorisPusherImpl_t>
struct MultiBorisFunctors
{
    static_assert(all_are<ParticleType>(particle_type));

    using GridLayout_t    = MultiBorisPusherImpl_t::GridLayout_t;
    using Particles_t     = MultiBorisPusherImpl_t::Particles_t;
    using Electromag_t    = MultiBorisPusherImpl_t::Electromag_t;
    using Interpolator_t  = MultiBorisPusherImpl_t::Interpolator_t;
    using ParticleArray_v = Particles_t::view_t;

    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Tile_vt       = Field_t::value_type::value_type;
    using VecField_vt   = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt = basic::Electromag<VecField_vt>;

    static constexpr auto dim = GridLayout_t::dimension;
    static_assert(Particles_t::storage_mode == StorageMode::VECTOR);

    MultiBorisFunctors(auto& in, auto& view, auto& pop, auto& parts, auto& em)
        : pps{*parts}
        , electromag{em}
        , dto2m{0.5 * in.dt / pop.mass()}
        , halfdt{in.mesh(view.layout.meshSize(), in.dt)}
    {
        check_particles(parts);
        check_particles_views(parts);
    }

    void check(auto const& particle) _PHARE_ALL_FN_
    {
        if constexpr (requires { check_particle(particle); })
            check_particle(particle);
    }

    void operator()(auto& in, [[maybe_unused]] auto const i)
    {
        if constexpr (Particles_t::alloc_mode == AllocatorMode::CPU)
            on_cpu_tiles(in);
        else if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            on_gpu_tiles(in, i);
        else
            throw std::runtime_error("Unknown alloc mode");
    }

    void on_gpu_tiles(auto& in, auto const i)
    {
#if !PHARE_HAVE_GPU
        throw std::runtime_error("NEEDS GPU IMPL!");
#else
        using Launcher = gpu::ChunkLauncher<false>;
        Launcher launcher{1, 0};
        launcher.b.x           = kernel::warp_size();
        launcher.g.x           = pps().size();
        auto const tile_picker = [pps = pps] __device__() {
            return std::make_tuple(blockIdx.x, &pps()[blockIdx.x], threadIdx.x,
                                   kernel::warp_size());
        };
        launcher.stream(in.streamer.streams[i],
                        [=, self = *this] __device__() mutable { self.per_tile(tile_picker); });
#endif
    }

    void on_cpu_tiles(auto& /*in*/)
    {
        std::size_t tileidx    = 0;
        auto const tile_picker = [&]() { return std::make_tuple(tileidx, &pps()[tileidx], 0, 1); };
        for (; tileidx < pps().size(); ++tileidx)
            per_tile(tile_picker);
    }

    void per_tile(auto const& tile_picker) _PHARE_ALL_FN_
    {
        auto&& [tile_idx, tileptr, tidx, ws] = tile_picker();
        auto& tile                           = *tileptr;
        auto const& layout                   = electromag.E[0][tile_idx].layout();
        auto const& em                       = em_tile(tile_idx);

        using enum LayoutMode;

        auto& parts          = tile();
        auto const tile_cell = pps.local_cell(tile.lower);

        auto constexpr static tracker = [](auto&&... args) _PHARE_ALL_FN_ {
            return make_particle_tracker<Particles_t::layout_mode, particle_type, dim>(args...);
        };

        if constexpr (Particles_t::layout_mode == AoSPCTS)
        {
            // GPU_UNIFIED for AoSPCTS to be handled separately later
            for (auto const& bix : parts.local_box())
            {
                auto& cell_particles = parts.particles_(bix);
                for (std::size_t pid = 0; pid < cell_particles.size(); ++pid)
                {
                    auto const& old_cell = cell_particles[pid].iCell();
                    auto const pt        = tracker(old_cell, pps.local_tile_cell(old_cell));
                    per_particle(cell_particles[pid], layout, pt, pid, em);
                }
            }
        }
        else
        {
            auto const each = pps()[tile_idx]().size() / ws;

            auto const one = [&](std::size_t const pidx) _PHARE_ALL_FN_ {
                per_any_particle(parts, layout, tracker(parts.iCell(pidx), tile_cell), pidx, em);
            };

            std::size_t pid = 0;
            for (; pid < each; ++pid)
                one(pid * ws + tidx);
            if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
                if (tidx < parts.size() - (ws * each))
                    one(pid * ws + tidx);
        }
    }

    void per_any_particle(auto& particles, auto&&... args) _PHARE_ALL_FN_
    {
        auto const& pidx = std::get<2>(std::forward_as_tuple(args...));
#if PHARE_HAVE_THRUST
        using enum LayoutMode;
        if constexpr (any_in(Particles_t::layout_mode, SoA, SoAPC, SoATS))
            per_particle(SoAZipParticle{particles, pidx}, args...);
        else
#endif
            per_particle(particles[pidx], args...);
    }


    void per_particle_still_in_ghost_box(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode             = Particles_t::alloc_mode;
        auto const& [particle, layout, pt, pidx, em] = std::forward_as_tuple(args...);

        check_electromag(em);
        check(particle);
        {
            Interpolator_t interp;
            boris::accelerate(particle, interp.m2p(particle, em, layout), dto2m);
        }
        check(particle);
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);
        check(particle);

        if constexpr (boris_mode == MultiBorisMode::REF)
        {
            using enum LayoutMode;
            if constexpr (Particles_t::layout_mode == AoSPCTS)
            {
                // AoSPCTS tracks per-cell buckets even in the tile ghost layer, so any
                // particle whose cell changed needs registering — domain and level ghost
                // alike. pt was built in per_tile, before advance() ran.
                pps.template move_check<particle_type>(pt, pidx, particle);
            }
            else if constexpr (particle_type == ParticleType::Domain)
            {
                if (isIn(particle, pps.box()))
                    pps.template move_check<particle_type>(pt, pidx, particle);
            }
        }

        check(particle);
    }

    void per_particle(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode             = Particles_t::alloc_mode;
        auto const& [particle, layout, pt, pidx, em] = std::forward_as_tuple(args...);

        check(particle);
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);
        check(particle);

        if constexpr (particle_type == ParticleType::Domain)
            per_particle_still_in_ghost_box(args...);
        else if constexpr (particle_type == ParticleType::LevelGhost)
        {
            if (isIn(particle, pps.ghost_box()))
                per_particle_still_in_ghost_box(args...);
            else
            {
                using enum LayoutMode;
                // left the ghost box on the first half-step — register the departure or
                // the bucket holding it goes stale; move_check resolves it as a deletion
                if constexpr (boris_mode == MultiBorisMode::REF
                              and Particles_t::layout_mode == AoSPCTS)
                    pps.template move_check<particle_type>(pt, pidx, particle);
            }
        }
        else
        {
            PHARE_ASSERT(false);
        }
    }


    void per_copy_of_cpu_tile_pcell(auto& boxings, auto& view, auto& pop, auto& parts)
    {
        auto const& patch_id      = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = parts.box();

        for (std::size_t tile_idx = 0; tile_idx < parts().size(); ++tile_idx)
        {
            auto& pctile         = parts()[tile_idx];
            auto& cell_particles = pctile();
            bool const is_border = patch_box * grow(pctile, 1) != grow(pctile, 1);
            auto const& layout   = electromag.E[0][tile_idx].layout();
            auto const tile_em   = em_tile(tile_idx);
            auto& rhoP           = pop.particleDensity()[tile_idx];
            auto& rhoC           = pop.chargeDensity()[tile_idx];
            auto F = pop.flux().template as<VecField_vt>([&](auto& c) { return c()[tile_idx](); });
            Interpolator_t interp;

            for (auto const& bix : cell_particles.local_box(cell_particles.box()))
            {
                for (auto p : cell_particles(bix))
                {
                    if constexpr (particle_type == ParticleType::LevelGhost)
                        if (!isIn(p, parts.ghost_box()))
                            continue;

                    p.iCell() = boris::advance<AllocatorMode::CPU>(p, halfdt);
                    boris::accelerate(p, interp.m2p(p, tile_em, layout), dto2m);
                    p.iCell() = boris::advance<AllocatorMode::CPU>(p, halfdt);

                    if constexpr (particle_type == ParticleType::Domain)
                    {
                        if (!is_border || isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                            interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                    }
                    else if constexpr (particle_type == ParticleType::LevelGhost)
                    {
                        if (isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                            interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                    }
                }
            }
        }
    }

    auto em_tile(auto const tidx) _PHARE_ALL_FN_
    {
        return electromag.template as<Electromag_vt>([&] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>([&](auto i) { return vf[i][tidx](); });
        });
    }

    void per_copy_of_cpu_tile(auto& boxings, auto& view, auto& pop, auto& parts)
    {
        auto const& patch_id = view.patchID();
        assert(boxings.count(patch_id));
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = parts.box();
        auto const is_border_     = [&](auto const& tile) {
            auto const tile_grow_box = grow(tile, 1);
            return patch_box * tile_grow_box != tile_grow_box;
        };

        for (std::size_t tile_idx = 0; tile_idx < parts().size(); ++tile_idx)
        {
            auto& tile           = parts()[tile_idx];
            auto const tile_cell = parts.local_cell(tile.lower);
            auto const& layout   = electromag.E[0][tile_idx].layout();
            auto const em        = em_tile(tile_idx);
            auto& rhoP           = pop.particleDensity()[tile_idx];
            auto& rhoC           = pop.chargeDensity()[tile_idx];
            auto F = pop.flux().template as<VecField_vt>([&](auto& c) { return c()[tile_idx](); });
            Interpolator_t interp;

            auto on_tile = [&]<bool border>() {
                for (auto p : tile())
                {
                    per_particle(p, layout, tile_cell, 0, em);
                    if constexpr (!border)
                        interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                    else
                    {
                        if (isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                            interp.particleToMesh(p, rhoP(), rhoC(), F, layout);
                    }
                }
            };

            if (is_border_(tile))
                on_tile.template operator()<true>();
            else
                on_tile.template operator()<false>();
        }
    }

    ParticleArray_v pps;
    Electromag_t::Super const electromag;
    double const dto2m;
    std::array<double, dim> halfdt;
};


// ── Common base: deduplicates move_rest / move_cpu_copy / sync_ref ───────────────────
// Derived must provide:
//   template<auto pt, auto mode> using Functors = MultiBorisFunctors<pt, mode, Derived>;
//   template<auto type> static void sync_particles(auto& particles, auto& stream);
//   using Particles_t = ...;

template<typename Derived>
struct MultiBorisPusherImplBase
{
    template<auto mode, typename ModelAccessor>
    static void move_rest(MultiBoris<ModelAccessor>& in, auto const i)
    {
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            domain.reset_views();
            typename Derived::template Functors<ParticleType::Domain, mode>{in, view, pop, domain,
                                                                            em}(in, i);

            auto& level_ghost = pop.levelGhostParticles();
            level_ghost.reset_views();
            typename Derived::template Functors<ParticleType::LevelGhost, mode>{
                in, view, pop, level_ghost, em}(in, i);
        }
    }

    template<auto mode, typename ModelAccessor>
    static void move_cpu_copy(MultiBoris<ModelAccessor>& in, auto& boxings, auto const i)
    {
        using enum LayoutMode;
        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        auto const per_parts = [&]<auto particle_type>(auto& pop, auto& parts) {
            parts.reset_views();
            typename Derived::template Functors<particle_type, mode> fns{in, view, pop, parts, em};
            if constexpr (any_in(Derived::Particles_t::layout_mode, AoSPCTS))
                fns.per_copy_of_cpu_tile_pcell(boxings, view, pop, parts);
            else
                fns.per_copy_of_cpu_tile(boxings, view, pop, parts);
        };

        for (auto& pop : ions)
        {
            per_parts.template operator()<ParticleType::Domain>(pop, pop.domainParticles());
            per_parts.template operator()<ParticleType::LevelGhost>(pop, pop.levelGhostParticles());
        }
    }

    template<typename ModelAccessor>
    static void sync_ref(MultiBoris<ModelAccessor>& in, auto const i)
    {
        if constexpr (Derived::Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            in.streamer.streams[i].sync();

        auto view      = in.accessor[i];
        auto [ions, _] = view.args;

        for (auto& pop : ions)
        {
            auto& domain = pop.domainParticles();
            Derived::template sync_particles<ParticleType::Domain>(domain, in.streamer.streams[i]);

            auto& level_ghost = pop.levelGhostParticles();
            Derived::template sync_particles<ParticleType::LevelGhost>(level_ghost,
                                                                       in.streamer.streams[i]);
        }
    }
};


// ── MultiBorisPusherImpl specialisations for AoSTS ────────────────────────────────────

template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImplBase<MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles,
                                                    Electromag, Interpolator>>
{
    using This  = MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                                       Interpolator>;
    using Super = MultiBorisPusherImplBase<This>;

public:
    static constexpr auto dim = GridLayout::dimension;
    using GridLayout_t        = GridLayout;
    using Particles_t         = Particles;
    using Electromag_t        = Electromag;
    using Interpolator_t      = Interpolator;


    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, This>;

    template<auto type>
    static void sync_particles(auto& particles, auto& stream)
    {
        sync_ts<type>(particles, stream);
    }


    template<MultiBorisMode mode = MultiBorisMode::REF, typename ModelAccessor>
    static void move(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        static constexpr auto copy   = mode == MultiBorisMode::COPY;
        static constexpr auto is_cpu = Particles_t::alloc_mode == AllocatorMode::CPU;
        static constexpr auto is_gpu = Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED;

        if constexpr (copy and is_gpu)
        {
            std::vector<default_span_size_t> sizes;
            sizes.reserve(in.accessor.size());
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
                sizes.push_back(boxings.at(in.accessor[i].patchID()).nonLevelGhostBox.size());
            in.gpu_nlgb = typename MultiBoris<ModelAccessor>::GpuBoxSpanSet_t{std::move(sizes)};
            std::size_t offset = 0;
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
            {
                auto const& nlgb = boxings.at(in.accessor[i].patchID()).nonLevelGhostBox;
                std::copy(nlgb.begin(), nlgb.end(), in.gpu_nlgb.vec.begin() + offset);
                offset += nlgb.size();
            }
        }

        auto move = [&](auto const i) mutable {
            if constexpr (copy and is_cpu)
                Super::template move_cpu_copy<mode>(in, boxings, i);
            else if constexpr (copy and is_gpu)
                move_gpu_copy<mode>(in, i);
            else
                Super::template move_rest<mode>(in, i);
        };
        auto sync = [&](auto const i) mutable {
            if constexpr (not copy)
                Super::sync_ref(in, i);
        };

        if constexpr (in.opts.use_main_thread)
        {
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
            {
                move(i);
                sync(i);
            }
        }
        else
        {
            // .host() takes a forwarding reference: passing the named lvalues directly
            // would deduce a reference type and store a dangling ref to this stack frame
            // once move() returns, so move them in to get a real, owned copy instead
            in.streamer.host(std::move(move));
            in.streamer.host(std::move(sync));
        }
    }


    template<auto mode, typename ModelAccessor>
    static void move_gpu_copy(MultiBoris<ModelAccessor>& in, auto const i)
    {
#if PHARE_HAVE_GPU
        using Tile_vt_ = Electromag_t::vecfield_type::field_type::value_type;
        using Interpolating_
            = Interpolating<Particles_t, GridLayout_t::interp_order,
                            (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)>;

        auto view       = in.accessor[i];
        auto [ions, em] = view.args;

        for (auto& pop : ions)
        {
            auto const dto2m      = 0.5 * in.dt / pop.mass();
            auto const halfdt     = MultiBoris<ModelAccessor>::mesh(view.layout.meshSize(), in.dt);
            auto const filter_box = in.gpu_nlgb[i];
            auto rhop             = pop.particleDensity();
            auto rhoc             = pop.chargeDensity();
            auto flux             = *pop.flux();
            auto const ds = static_cast<std::uint32_t>(ions.chargeDensity().max_tile_size());

            auto const launch = [&](auto parts) {
                if (parts().size() == 0)
                    return;
                using Launcher = gpu::ChunkLauncher<false>;
                Launcher launcher{1, 0};
                launcher.b.x = kernel::warp_size();
                launcher.g.x = parts().size();
                launcher.ds
                    = ds * 5 * 8
                      + 2 * static_cast<std::uint32_t>(kernel::warp_size())
                            * static_cast<std::uint32_t>(sizeof(typename Particles_t::Particle_t))
                      + static_cast<std::uint32_t>(sizeof(int));
                assert(launcher.ds < 65000);
                launcher.stream(in.streamer.streams[i], [=] __device__() mutable {
                    Interpolating_::template on_tiles_push_deposit<Tile_vt_, ParticleType::Domain,
                                                                   Particles_t::alloc_mode>(
                        parts, em, flux, rhop, rhoc, filter_box, dto2m, halfdt);
                });
            };

            launch(*pop.domainParticles());
            launch(*pop.levelGhostParticles());
        }
#endif
    }
};


// AoSMapped: reuse AoSTS tile-based functors
template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSMapped, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImpl<LayoutMode::AoSTS, alloc, GridLayout, Particles, Electromag,
                           Interpolator>
{
};


// AoSPCTS: tiles of per-cell AoS particles, tiled fields (same as AoSTS)
// Particle iteration differs: each tile contains a PerCellParticles container
// with flat iterators over all particles in that tile.
template<AllocatorMode alloc, typename GridLayout, typename Particles, typename Electromag,
         typename Interpolator>
struct MultiBorisPusherImpl<LayoutMode::AoSPCTS, alloc, GridLayout, Particles, Electromag,
                            Interpolator>
    : MultiBorisPusherImplBase<MultiBorisPusherImpl<LayoutMode::AoSPCTS, alloc, GridLayout,
                                                    Particles, Electromag, Interpolator>>
{
    using This = MultiBorisPusherImpl<LayoutMode::AoSPCTS, alloc, GridLayout, Particles, Electromag,
                                      Interpolator>;
    using Super = MultiBorisPusherImplBase<This>;

    static constexpr auto dim = GridLayout::dimension;


    using GridLayout_t   = GridLayout;
    using Particles_t    = Particles;
    using Electromag_t   = Electromag;
    using Interpolator_t = Interpolator;

    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Tile_vt       = Field_t::value_type::value_type;
    using VecField_vt   = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt = basic::Electromag<VecField_vt>;

    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, This>;

    template<auto type>
    static void sync_particles(auto& particles, auto& stream)
    {
        sync_pc_ts<type>(particles, stream);
    }

    template<MultiBorisMode mode = MultiBorisMode::REF, typename ModelAccessor>
    static void move(MultiBoris<ModelAccessor>& in, auto const& boxings)
    {
        static constexpr auto copy   = mode == MultiBorisMode::COPY;
        static constexpr auto is_cpu = Particles_t::alloc_mode == AllocatorMode::CPU;

        auto move = [&](auto const i) mutable {
            if constexpr (copy and is_cpu)
                Super::template move_cpu_copy<mode>(in, boxings, i);
            else
                Super::template move_rest<mode>(in, i);
        };
        auto sync = [&](auto const i) mutable {
            if constexpr (not copy)
                Super::sync_ref(in, i);
        };

        if constexpr (in.opts.use_main_thread)
        {
            for (std::size_t i = 0; i < in.accessor.size(); ++i)
            {
                move(i);
                sync(i);
            }
        }
        else
        {
            // see AoSTS::move() above: .host() must receive an rvalue or it stores a
            // dangling reference to this stack frame instead of an owned copy
            in.streamer.host(std::move(move));
            in.streamer.host(std::move(sync));
        }
    }
};

} // namespace PHARE::core


#endif
