#ifndef PHARE_CORE_PUSHER_BORIS_DETAIL_MULTI_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_DETAIL_MULTI_BORIS_HPP

#include "core/utilities/span.hpp"
#include "core/utilities/kernels.hpp"
#include "core/utilities/thread_pool.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

namespace PHARE::core
{

enum class MultiBorisMode : std::uint16_t { REF = 0, COPY };

struct MultiBorisOptions
{
    bool use_main_thread = false; // true for perf
};

template<auto particle_type, auto boris_mode, typename Backend_t>
struct MultiBorisFunctors;

template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
struct MultiBorisBackend;


// ── MultiBoris state struct / ModelAccessor ───────────────────────────────────────────
// default dispatch: BS thread pool (or sequential with MultiBorisOptions::use_main_thread)
// see core::mkn_xyz::MultiBoris below for the mkn.gpu ThreadedStreamLauncher dispatch

template<typename ModelAccessor, typename Interpolator>
struct MultiBoris
{
    static constexpr auto dim = ModelAccessor::GridLayout_t::dimension;
    using Model_t             = ModelAccessor::Model_t;
    using GridLayout_t        = ModelAccessor::GridLayout_t;
    using ParticleArray_t     = Model_t::particle_array_type;
    using Electromag_t        = Model_t::electromag_type;
    using Vecfield_t          = Electromag_t::vecfield_type;
    using Field_t             = Vecfield_t::field_type;
    using ParticleArray_v     = ParticleArray_t::view_t;
    using Backend = MultiBorisBackend<ParticleArray_t::layout_mode, ModelAccessor, Interpolator>;

    MultiBoris(double const dt_, ModelAccessor& _accessor)
        : dt{dt_}
        , accessor{_accessor}
    {
    }

    template<MultiBorisMode mode = MultiBorisMode::REF>
    void move(auto const& boxings);

    // No thread pools at all: patches then tiles, both fully sequential on the caller's
    // thread. Used when MultiBorisOptions::use_main_thread is set (perf comparisons /
    // environments where spinning up pools isn't worth it).
    template<MultiBorisMode mode>
    void move_sequential(auto const& boxings);

    // General case: two levels of parallelism. Each patch is assigned one pool
    // (round-robin over whichever is next free), and within that pool every tile of
    // the patch is its own task -- so with N pools of M threads, up to N patches and
    // N*M tiles are in flight simultaneously. Dispatch happens entirely from this
    // (calling) thread so no pool ever detach_task-then-waits on itself.
    template<MultiBorisMode mode>
    void move_pooled(auto const& boxings);

    double const dt;
    ModelAccessor& accessor;

    auto static mesh(std::array<double, dim> const& ms, double const& ts)
    {
        std::array<double, dim> halfDtOverDl;
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl),
                       [ts](double const& x) { return 0.5 * ts / x; });
        return halfDtOverDl;
    }
};


template<typename ModelAccessor, typename Interpolator>
template<MultiBorisMode mode>
void MultiBoris<ModelAccessor, Interpolator>::move(auto const& boxings)
{
    if constexpr (MultiBorisOptions{}.use_main_thread)
        move_sequential<mode>(boxings);
    else
        move_pooled<mode>(boxings);
}


template<typename ModelAccessor, typename Interpolator>
template<MultiBorisMode mode>
void MultiBoris<ModelAccessor, Interpolator>::move_sequential(auto const& boxings)
{
    static constexpr auto copy   = mode == MultiBorisMode::COPY;
    static constexpr auto is_cpu = ParticleArray_t::alloc_mode == AllocatorMode::CPU;

    for (std::size_t i = 0; i < accessor.size(); ++i)
    {
        if constexpr (copy and is_cpu)
            Backend::template move_cpu_copy<mode>(*this, boxings, i);
        else
            Backend::template move_rest<mode>(*this, i);

        if constexpr (not copy)
            Backend::sync_ref(*this, i);
    }
}


template<typename ModelAccessor, typename Interpolator>
template<MultiBorisMode mode>
void MultiBoris<ModelAccessor, Interpolator>::move_pooled(auto const& boxings)
{
    static constexpr auto copy   = mode == MultiBorisMode::COPY;
    static constexpr auto is_cpu = ParticleArray_t::alloc_mode == AllocatorMode::CPU;

    auto& TP = ThreadPool::INSTANCE();

    for (std::size_t i = 0; i < accessor.size(); ++i)
    {
        auto& pool = TP.get_pool(TP.first_ready_idx());
        if constexpr (copy and is_cpu)
            Backend::template move_cpu_copy_pooled<mode>(pool, *this, boxings, i);
        else
            Backend::template move_rest_pooled<mode>(pool, *this, i);
    }
    TP.sync(); // every patch's tiles, across every pool, are done

    if constexpr (not copy)
    {
        for (std::size_t i = 0; i < accessor.size(); ++i)
        {
            auto& pool = TP.get_pool(TP.first_ready_idx());
            pool.detach_task([this, i] { Backend::sync_ref(*this, i); });
        }
        TP.sync();
    }
}


// ── Per patch steps, shared by every dispatch (core::MultiBoris / mkn_xyz::MultiBoris) ────
// tiles of particles, tiled fields: AoSPCTS (per-cell AoS per tile), AoSTS/AoSMapped (flat)
// `in` is either MultiBoris state struct, GPU paths require mkn_xyz::MultiBoris (in.streamer)

template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
struct MultiBorisBackend
{
    using MultiBoris_t   = MultiBoris<ModelAccessor, Interpolator>;
    using GridLayout_t   = MultiBoris_t::GridLayout_t;
    using Particles_t    = MultiBoris_t::ParticleArray_t;
    using Electromag_t   = MultiBoris_t::Electromag_t;
    using Interpolator_t = Interpolator;

    static constexpr auto dim = GridLayout_t::dimension;

    // only flat tiles have a push+deposit GPU kernel, AoSPCTS GPU COPY uses move_rest
    static constexpr bool has_gpu_copy = any_in(layout, LayoutMode::AoSTS, LayoutMode::AoSMapped);

    template<auto pt, auto mode>
    using Functors = MultiBorisFunctors<pt, mode, MultiBorisBackend>;

    template<auto type>
    static void sync_particles(auto& particles, auto&... stream)
    {
        particles.template on_moved<type>(stream...);
    }

    template<auto mode>
    static void move_rest(auto& in, auto const i);

    template<auto mode>
    static void move_cpu_copy(auto& in, auto& boxings, auto const i);

    // ── Thread-pooled variants: level 1 (patch -> pool) is dispatched by the caller
    // (move_pooled), which hands us the specific pool already assigned to patch `i`.
    // Here we do level 2: one task per tile, submitted to that SAME pool so all of its
    // threads_per_pool threads work this patch's tiles at once. These must be called
    // from the dispatching thread itself, never from inside a task already running on
    // `pool` — detach_task-then-wait on one's own pool deadlocks.

    template<auto mode>
    static void move_rest_pooled(auto& pool, auto& in, auto const i);

    template<auto mode>
    static void move_cpu_copy_pooled(auto& pool, auto& in, auto& boxings, auto const i);

    static void sync_ref(auto& in, auto const i);

    // GPU only: copies every patch's nonLevelGhostBox into managed memory for the kernels
    static void prepare_gpu_copy(auto& in, auto const& boxings);

    template<auto mode>
    static void move_gpu_copy(auto& in, auto const i);

    // Combined Boris push + particle-to-mesh deposit in a single GPU kernel.
    // Originals are untouched: each thread works on a register copy.
    // Qualifying particles (isIn filter_box and not in exclude_boxes after push) are atomically
    // placed into a shared warp buffer. Once the buffer holds a full warp the whole warp
    // deposits together; any remainder is flushed at the end.
    //
    // Shared memory layout (launcher.ds must cover):
    //   [4 * ziz * sizeof(double)]            field accumulators
    //   [2 * warp_size * sizeof(Particle_t)]  warp buffer
    //   [sizeof(int)]                          atomic fill counter
    //
    // 2*warp_size slots are the minimum needed: before a batch buf_count can be
    // warp_size-1, and if every thread in the batch qualifies buf_count reaches
    // 2*warp_size-1 — still within bounds.
    template<typename Field_t, auto particle_type_v, auto alloc_mode_v>
    static void on_tiles_push_deposit(auto& particles, auto const& electromag, auto& Flux,
                                      auto& Rhop, auto const& filter_box, auto const& exclude_boxes,
                                      double const coef1, auto const& halfdt,
                                      double const coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        static_assert(
            std::is_same_v<Interpolator_t,
                           core::Interpolator<dim, GridLayout_t::options.interp_order, true>>,
            "GPU must be atomic");

        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;
        using Scalar_vt     = typename Field_t::value_type;
        using VecField_vt   = basic::TensorField<Scalar_vt, 1>;
        using Electromag_vt = basic::Electromag<VecField_vt>;
        using Particle_t    = ParticleDefaults<dim>::Particle_t;

        extern __shared__ double data[];

        auto& rhop_tile = Rhop[blockIdx.x];
        auto ft         = Flux.template as<TensorField_t>([&](auto& c) { return c[blockIdx.x](); });
        auto& pdensity  = rhop_tile();

        auto const tbox                   = rhop_tile.ghost_box();
        auto const ziz                    = product(tbox.shape());
        auto const tidx                   = threadIdx.x;
        auto const ws                     = kernel::warp_size();
        auto const& [flux0, flux1, flux2] = ft();
        auto const r0                     = pdensity;
        auto flux                         = ft;

        auto rhop = make_array_view(&data[ziz * 0], pdensity.shape());
        auto fx   = make_array_view(&data[ziz * 1], flux[0].shape());
        auto fy   = make_array_view(&data[ziz * 2], flux[1].shape());
        auto fz   = make_array_view(&data[ziz * 3], flux[2].shape());

        pdensity.reset(rhop);
        flux[0].reset(fx);
        flux[1].reset(fy);
        flux[2].reset(fz);

        // Warp buffer and atomic counter after the field region.
        auto* const part_buf = reinterpret_cast<Particle_t*>(reinterpret_cast<char*>(&data[0])
                                                             + ziz * 4 * sizeof(double));
        auto& buf_count      = *reinterpret_cast<int*>(reinterpret_cast<char*>(part_buf)
                                                       + 2 * ws * sizeof(Particle_t));

        for (int i = tidx; i < static_cast<int>(ziz * 4); i += blockDim.x)
            data[i] = 0;
        if (tidx == 0)
            buf_count = 0;

        __syncthreads();

        // Per-tile EM view and layout
        auto const em = electromag.template as<Electromag_vt>([=] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>(
                [&vf] _PHARE_ALL_FN_(auto i) { return vf[i][blockIdx.x](); });
        });
        auto const& tile_layout = electromag.E[0][blockIdx.x].layout();

        auto& tile        = particles()[blockIdx.x];
        auto& parts       = tile();
        auto const n      = parts.size();
        auto const nbatch = (n + ws - 1) / ws;

        for (std::size_t batch = 0; batch < nbatch; ++batch)
        {
            std::size_t const pidx = batch * ws + tidx;
            if (pidx < n)
            {
                Particle_t p = parts[pidx]; // register copy — original untouched

                p.iCell() = boris::advance<alloc_mode_v>(p, halfdt);
                {
                    Interpolator_t interp;
                    boris::accelerate(p, interp.m2p(p, em, tile_layout), coef1);
                }
                p.iCell() = boris::advance<alloc_mode_v>(p, halfdt);

                if (isIn(p.iCell(), filter_box) and !isIn(p.iCell(), exclude_boxes))
                {
                    int const slot = atomicAdd(&buf_count, 1);
                    part_buf[slot] = p;
                }
            }

            __syncthreads();

            // Once a full warp is available, deposit it and compact the buffer.
            if (buf_count >= static_cast<int>(ws))
            {
                Interpolator_t interp;
                interp.particleToMesh(part_buf[tidx], pdensity, flux, tbox, coef);

                __syncthreads();

                int const remaining = buf_count - static_cast<int>(ws);
                if (tidx < remaining)
                    part_buf[tidx] = part_buf[tidx + ws];
                if (tidx == 0)
                    buf_count = remaining;

                __syncthreads();
            }
        }

        // Flush any remaining particles (partial warp).
        if (tidx < buf_count)
        {
            Interpolator_t interp;
            interp.particleToMesh(part_buf[tidx], pdensity, flux, tbox, coef);
        }

        __syncthreads();

        pdensity.reset(r0);

        for (int i = tidx; i < static_cast<int>(ziz); i += blockDim.x)
        {
            pdensity.data()[i] += rhop.data()[i];
            flux0.data()[i] += fx.data()[i];
            flux1.data()[i] += fy.data()[i];
            flux2.data()[i] += fz.data()[i];
        }
#endif
    }
};


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
template<auto mode>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::move_rest(auto& in, auto const i)
{
    auto view       = in.accessor[i];
    auto [ions, em] = view.args;

    for (auto& pop : ions)
    {
        auto& domain = pop.domainParticles();
        domain.reset_views();
        (Functors<ParticleType::Domain, mode>{in, view, pop, domain, em})(in, i);

        auto& level_ghost = pop.levelGhostParticles();
        level_ghost.reset_views();
        (Functors<ParticleType::LevelGhost, mode>{in, view, pop, level_ghost, em})(in, i);
    }
}


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
template<auto mode>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::move_cpu_copy(auto& in, auto& boxings,
                                                                           auto const i)
{
    auto view       = in.accessor[i];
    auto [ions, em] = view.args;

    auto const per_parts = [&]<auto particle_type>(auto& pop, auto& parts) {
        parts.reset_views();
        Functors<particle_type, mode> fns{in, view, pop, parts, em};
        fns.per_copy_of_cpu_tile(boxings, view, pop);
    };

    for (auto& pop : ions)
    {
        per_parts.template operator()<ParticleType::Domain>(pop, pop.domainParticles());
        per_parts.template operator()<ParticleType::LevelGhost>(pop, pop.levelGhostParticles());
    }
}


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
template<auto mode>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::move_rest_pooled(auto& pool, auto& in,
                                                                              auto const i)
{
    auto view       = in.accessor[i];
    auto [ions, em] = view.args;

    for (auto& pop : ions)
    {
        auto& domain = pop.domainParticles();
        domain.reset_views();
        auto& level_ghost = pop.levelGhostParticles();
        level_ghost.reset_views();

        using DomainFn = Functors<ParticleType::Domain, mode>;
        using GhostFn  = Functors<ParticleType::LevelGhost, mode>;

        // shared_ptr keeps the (self-contained: own pps view + electromag copy)
        // Functors alive for every deferred per-tile task below.
        auto domain_fn = std::make_shared<DomainFn>(in, view, pop, domain, em);
        auto ghost_fn  = std::make_shared<GhostFn>(in, view, pop, level_ghost, em);

        auto const n_tiles = domain_fn->pps().size();
        assert(n_tiles == ghost_fn->pps().size());

        for (std::size_t t = 0; t < n_tiles; ++t)
            pool.detach_task([domain_fn, ghost_fn, t] {
                domain_fn->one_tile(t);
                ghost_fn->one_tile(t);
            });
    }
}


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
template<auto mode>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::move_cpu_copy_pooled(auto& pool,
                                                                                  auto& in,
                                                                                  auto& boxings,
                                                                                  auto const i)
{
    using DomainFn = Functors<ParticleType::Domain, mode>;
    using GhostFn  = Functors<ParticleType::LevelGhost, mode>;

    auto view       = in.accessor[i];
    auto [ions, em] = view.args;

    for (auto& pop : ions)
    {
        auto& domain = pop.domainParticles();
        domain.reset_views();
        auto& level_ghost = pop.levelGhostParticles();
        level_ghost.reset_views();

        auto domain_fn     = std::make_shared<DomainFn>(in, view, pop, domain, em);
        auto ghost_fn      = std::make_shared<GhostFn>(in, view, pop, level_ghost, em);
        auto const n_tiles = domain_fn->pps().size();
        assert(n_tiles == ghost_fn->pps().size());

        // view/pop are cheap handles, safe to capture by copy. Domain and LevelGhost
        // for the SAME tile must run in the SAME task (shared rhoP/flux);
        // different tiles are independent.
        for (std::size_t t = 0; t < n_tiles; ++t)
            pool.detach_task([domain_fn, ghost_fn, t, view, pop, &boxings]() mutable {
                domain_fn->one_copy_tile(t, boxings, view, pop);
                ghost_fn->one_copy_tile(t, boxings, view, pop);
            });
    }
}


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::sync_ref(auto& in, auto const i)
{
    static constexpr bool has_streams = requires { in.streamer; };

    if constexpr (has_streams and Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
        in.streamer.streams[i].sync();

    auto view      = in.accessor[i];
    auto [ions, _] = view.args;

    for (auto& pop : ions)
    {
        auto& domain      = pop.domainParticles();
        auto& level_ghost = pop.levelGhostParticles();

        if constexpr (has_streams)
        {
            sync_particles<ParticleType::Domain>(domain, in.streamer.streams[i]);
            sync_particles<ParticleType::LevelGhost>(level_ghost, in.streamer.streams[i]);
        }
        else
        {
            sync_particles<ParticleType::Domain>(domain);
            sync_particles<ParticleType::LevelGhost>(level_ghost);
        }
    }
}


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::prepare_gpu_copy(auto& in,
                                                                              auto const& boxings)
{
    using GpuBoxSpanSet_t = std::remove_cvref_t<decltype(in.gpu_nlgb)>;

    auto const to_gpu = [&](auto const& patch_boxes) {
        std::vector<default_span_size_t> sizes;
        sizes.reserve(in.accessor.size());
        for (std::size_t i = 0; i < in.accessor.size(); ++i)
            sizes.push_back(patch_boxes(boxings.at(in.accessor[i].patchID())).size());
        GpuBoxSpanSet_t gpu_boxes{std::move(sizes)};
        std::size_t offset = 0;
        for (std::size_t i = 0; i < in.accessor.size(); ++i)
        {
            auto const& boxes = patch_boxes(boxings.at(in.accessor[i].patchID()));
            std::copy(boxes.begin(), boxes.end(), gpu_boxes.vec.begin() + offset);
            offset += boxes.size();
        }
        return gpu_boxes;
    };

    in.gpu_nlgb      = to_gpu([](auto const& b) -> auto& { return b.nonLevelGhostBox; });
    in.gpu_domain    = to_gpu([](auto const& b) { return std::vector{b.domainBox}; });
    in.gpu_non_owned = to_gpu([](auto const& b) -> auto& { return b.nonOwnedBoxes; });
}


template<LayoutMode layout, typename ModelAccessor, typename Interpolator>
template<auto mode>
void MultiBorisBackend<layout, ModelAccessor, Interpolator>::move_gpu_copy(
    [[maybe_unused]] auto& in, [[maybe_unused]] auto const i)
{
#if PHARE_HAVE_GPU
    static_assert(has_gpu_copy);
    using Tile_vt_ = Electromag_t::vecfield_type::field_type::value_type;

    auto view       = in.accessor[i];
    auto [ions, em] = view.args;

    for (auto& pop : ions)
    {
        auto const coef1      = pop.charge() * (0.5 * in.dt / pop.mass());
        auto const halfdt     = in.mesh(view.layout.meshSize(), in.dt);
        auto const nlgb      = in.gpu_nlgb[i];
        auto const domain    = in.gpu_domain[i];
        auto const non_owned = in.gpu_non_owned[i];
        auto rhop             = pop.particleDensity();
        auto flux             = *pop.flux();
        auto const ds         = static_cast<std::uint32_t>(ions.chargeDensity().max_tile_size());

        auto const launch = [&](auto parts, auto const filter_box, auto const exclude_boxes) {
            if (parts().size() == 0)
                return;
            using Launcher = gpu::ChunkLauncher<false>;
            Launcher launcher{1, 0};
            launcher.b.x = kernel::warp_size();
            launcher.g.x = parts().size();
            launcher.ds
                = ds * 4 * 8
                  + 2 * static_cast<std::uint32_t>(kernel::warp_size())
                        * static_cast<std::uint32_t>(sizeof(typename Particles_t::Particle_t))
                  + static_cast<std::uint32_t>(sizeof(int));
            assert(launcher.ds < 65000);
            launcher.stream(in.streamer.streams[i], [=] __device__() mutable {
                on_tiles_push_deposit<Tile_vt_, ParticleType::Domain, Particles_t::alloc_mode>(
                    parts, em, flux, rhop, filter_box, exclude_boxes, coef1, halfdt);
            });
        };

        launch(*pop.domainParticles(), nlgb, decltype(non_owned){});
        launch(*pop.levelGhostParticles(), domain, non_owned);
    }
#else
    throw std::runtime_error("NEEDS GPU IMPL!");
#endif
}



// ── Per-tile functors (shared between CPU and GPU, and every dispatch) ─────────────────
// per-tile containers differ per layout (per-cell for AoSPCTS, flat cell-mapped for
// AoSCMTS, flat for AoSTS): each is handled by its own requires(layout) overload

template<auto particle_type, auto boris_mode, typename Backend_t>
struct MultiBorisFunctors
{
    static_assert(all_are<ParticleType>(particle_type));

    using GridLayout_t    = Backend_t::GridLayout_t;
    using Particles_t     = Backend_t::Particles_t;
    using Electromag_t    = Backend_t::Electromag_t;
    using Interpolator_t  = Backend_t::Interpolator_t;
    using ParticleArray_v = Particles_t::view_t;

    using Vecfield_t    = Electromag_t::vecfield_type;
    using Field_t       = Vecfield_t::field_type;
    using Tile_vt       = Field_t::value_type::value_type;
    using VecField_vt   = basic::TensorField<Tile_vt, 1>;
    using Electromag_vt = basic::Electromag<VecField_vt>;

    static constexpr auto dim = GridLayout_t::dimension;
    using Particle_t          = ParticleDefaults<dim>::Particle_t;

    static constexpr std::size_t copy_batch_size = 1024;
    static_assert(Particles_t::storage_mode == StorageMode::VECTOR);

    MultiBorisFunctors(auto& in, auto& view, auto& pop, auto& parts, auto& em)
        : pps{*parts}
        , electromag{em}
        , coef1{pop.charge() * (0.5 * in.dt / pop.mass())}
        , halfdt{in.mesh(view.layout.meshSize(), in.dt)}
        , particles{parts}
    {
        check_particles(parts);
        check_particles_views(parts);
    }

    void operator()(auto& in, [[maybe_unused]] auto const i)
    {
        // GPU_UNIFIED without streams (core::MultiBoris) falls back to the CPU loop
        if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED
                      and requires { in.streamer; })
            on_gpu_tiles(in, i);
        else
            on_cpu_tiles();
    }

    void on_gpu_tiles([[maybe_unused]] auto& in, [[maybe_unused]] auto const i)
    {
#if !PHARE_HAVE_GPU
        throw std::runtime_error("NEEDS GPU IMPL!");
#else
        if constexpr (Particles_t::layout_mode == LayoutMode::AoSPCTS)
            on_gpu_tile_cells(in, i);
        else
            on_gpu_flat_tiles(in, i);
#endif
    }

    // AoSPCTS: block = (tile, cell), block threads stride the cell's particles
    void on_gpu_tile_cells([[maybe_unused]] auto& in, [[maybe_unused]] auto const i)
    {
#if PHARE_HAVE_GPU
        using Launcher = gpu::TileCellLauncher<false>;

        auto const max_cells = Launcher::max_cells(pps);
        if (max_cells == 0)
            return;

        std::size_t const threads = kernel::warp_size();
        Launcher{pps().size(), max_cells, threads}.stream(
            in.streamer.streams[i], [self = *this, threads] __device__() mutable {
                auto const tile_idx = Launcher::tile_idx();
                auto& tile          = self.pps()[tile_idx];
                auto& parts         = tile();
                auto const& lobox   = parts.local_box();
                if (Launcher::cell_idx() >= lobox.size())
                    return;
                auto const& bix = *(lobox.begin() + Launcher::cell_idx());
                self.per_cell(parts.particles_(bix), self.electromag.E[0][tile_idx].layout(),
                              self.em_tile(tile_idx), self.pps.local_cell(tile.lower),
                              Launcher::thread_idx(), threads);
            });
#endif
    }

    // flat tiles: block = tile, warp of threads striding the tile's particles
    void on_gpu_flat_tiles([[maybe_unused]] auto& in, [[maybe_unused]] auto const i)
    {
#if PHARE_HAVE_GPU
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

    void on_cpu_tiles()
    {
        for (std::size_t tileidx = 0; tileidx < pps().size(); ++tileidx)
            one_tile(tileidx);
    }

    // one tile's worth of work — the unit of dispatch for pool.detach_task() in the
    // thread-pooled path (see MultiBorisBackend::move_rest_pooled).
    void one_tile(std::size_t const tile_idx)
    {
        auto const tile_picker
            = [&]() { return std::make_tuple(tile_idx, &pps()[tile_idx], 0, 1); };
        per_tile(tile_picker);
    }

    static auto tracker(auto&&... args) _PHARE_ALL_FN_
    {
        return make_particle_tracker<Particles_t::layout_mode, particle_type, dim>(args...);
    }

    // tile_cell names the tile PHYSICALLY holding a particle (for level ghosts this is
    // also their cell's clamp owner)

    void per_tile(auto const& tile_picker)
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    _PHARE_ALL_FN_
    {
        auto&& [tile_idx, tileptr, tidx, ws] = tile_picker();
        auto& tile                           = *tileptr;
        auto const& layout                   = electromag.E[0][tile_idx].layout();
        auto const& em                       = em_tile(tile_idx);
        auto& parts                          = tile();
        auto const tile_cell                 = pps.local_cell(tile.lower);

        // GPU uses on_gpu_tile_cells, one block per cell
        for (auto const& bix : parts.local_box())
            per_cell(parts.particles_(bix), layout, em, tile_cell, tidx, ws);
    }

    // flat tiles (AoSCMTS, AoSTS, ...): pid is the index in the tile's array
    void per_tile(auto const& tile_picker)
        requires(Particles_t::layout_mode != LayoutMode::AoSPCTS)
    _PHARE_ALL_FN_
    {
        auto&& [tile_idx, tileptr, tidx, ws] = tile_picker();
        auto& tile                           = *tileptr;
        auto const& layout                   = electromag.E[0][tile_idx].layout();
        auto const& em                       = em_tile(tile_idx);
        auto& parts                          = tile();
        auto const tile_cell                 = pps.local_cell(tile.lower);

        auto const each = pps()[tile_idx]().size() / ws;

        auto const one = [&] _PHARE_ALL_FN_(std::size_t const pidx) {
            per_any_particle(parts, layout, tracker(parts.iCell(pidx), tile_cell), pidx, em);
        };

        std::size_t pid = 0;
        for (; pid < each; ++pid)
            one(pid * ws + tidx);
        if constexpr (Particles_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
            if (tidx < parts.size() - (ws * each))
                one(pid * ws + tidx);
    }


    // per-cell buckets are tracked even in the tile ghost layer, so any particle whose
    // cell changed needs registering — domain and level ghost alike
    void move_check(auto const& pt, std::size_t const pidx, auto& particle)
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    _PHARE_ALL_FN_
    {
        pps.template move_check<particle_type>(pt, pidx, particle);
    }

    // AoSCMTS registers on the vector: a cellmap can't be updated through its view
    void move_check(auto const& pt, std::size_t const pidx, auto& particle)
        requires(Particles_t::layout_mode == LayoutMode::AoSCMTS)
    {
        particles.template move_check<particle_type>(pt, pidx, particle);
    }

    // other flat tiles: only domain particles staying in the patch box register
    void move_check(auto const& pt, std::size_t const pidx, auto& particle)
        requires(not any_in(Particles_t::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
    _PHARE_ALL_FN_
    {
        if constexpr (particle_type == ParticleType::Domain)
            if (isIn(particle, pps.box()))
                pps.template move_check<particle_type>(pt, pidx, particle);
    }

    // AoSPCTS: one cell of a tile, particles from tidx strided by ws (CPU: 0, 1)
    void per_cell(auto& cell_particles, auto const& layout, auto const& em, auto const& tile_cell,
                  std::size_t const tidx, std::size_t const ws) _PHARE_ALL_FN_
    {
        for (std::size_t pid = tidx; pid < cell_particles.size(); pid += ws)
            per_particle(cell_particles[pid], layout,
                         tracker(cell_particles[pid].iCell(), tile_cell), pid, em);
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
        {
            Interpolator_t interp;
            boris::accelerate(particle, interp.m2p(particle, em, layout), coef1);
        }
        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);

        if constexpr (boris_mode == MultiBorisMode::REF)
            move_check(pt, pidx, particle); // pt was built in per_tile, before advance() ran
    }

    void per_particle(auto&&... args) _PHARE_ALL_FN_
    {
        static constexpr auto alloc_mode             = Particles_t::alloc_mode;
        auto const& [particle, layout, pt, pidx, em] = std::forward_as_tuple(args...);

        particle.iCell() = boris::advance<alloc_mode>(particle, halfdt);

        if constexpr (particle_type == ParticleType::Domain)
            per_particle_still_in_ghost_box(args...);
        else if constexpr (particle_type == ParticleType::LevelGhost)
        {
            if (isIn(particle, pps.ghost_box()))
                per_particle_still_in_ghost_box(args...);
            else if constexpr (boris_mode == MultiBorisMode::REF)
                // left the ghost box on the first half-step — register the departure or
                // the bucket holding it goes stale; move_check resolves it as a deletion
                move_check(pt, pidx, particle);
        }
        else
        {
            PHARE_ASSERT(false);
        }
    }

    // one tile's worth of the COPY-mode deposit pass. Domain and LevelGhost write into
    // the SAME tile's rhoP/flux, so callers must run this tile's Domain and
    // LevelGhost passes in one task (never two concurrent tasks for the same tile_idx)
    // — see MultiBorisBackend::move_cpu_copy_pooled. Different tiles are independent.
    void one_copy_tile(std::size_t const tile_idx, auto& boxings, auto& view, auto& pop)
        requires(any_in(Particles_t::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
    {
        auto const& patch_id      = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = pps.box();

        auto& tile           = pps()[tile_idx];
        bool const is_border = patch_box * grow(tile, 1) != grow(tile, 1);
        auto const& layout   = electromag.E[0][tile_idx].layout();
        auto const tile_em   = em_tile(tile_idx);
        auto& rhoP           = pop.particleDensity()[tile_idx];
        auto F               = tile_at(pop.flux(), tile_idx);
        Interpolator_t interp;

        auto const deposits = [&](auto const& p) {
            if constexpr (particle_type == ParticleType::Domain)
                return !is_border || isIn(p.iCell(), patch_boxings.nonLevelGhostBox);
            else
                // level ghosts live only in their clamp-owner tile, so each is seen
                // once: deposit any entering the patch domain into this tile's fields,
                // the tile reduction folds halo contributions into the neighbours
                return patch_boxings.isOwned(p.iCell());
        };

        // batched in phases (copy + half step, m2p + accelerate + half step, deposit) as
        // interleaving m2p and p2m per particle is slower, see updater_realloc's
        // IonUpdater1::Pusher::move_domain_copy
        std::array<Particle_t, copy_batch_size> batch;
        std::size_t n = 0;

        auto const flush = [&]() {
            std::size_t m = 0;
            for (std::size_t j = 0; j < n; ++j)
            {
                auto& p = batch[j];
                boris::accelerate(p, interp.m2p(p, tile_em, layout), coef1);
                p.iCell() = boris::advance<AllocatorMode::CPU>(p, halfdt);
                if (deposits(p))
                    batch[m++] = p;
            }
            for (std::size_t j = 0; j < m; ++j)
                interp.particleToMesh(batch[j], rhoP(), F, layout);
            n = 0;
        };

        auto const per_array = [&](auto const& ps) {
            for (auto const& p : ps)
            {
                if constexpr (particle_type == ParticleType::LevelGhost)
                    if (!isIn(p, pps.ghost_box()))
                        continue;

                batch[n]         = p;
                batch[n].iCell() = boris::advance<AllocatorMode::CPU>(batch[n], halfdt);
                if (++n == copy_batch_size)
                    flush();
            }
        };

        for_each_copy_array(tile(), per_array);
        flush();
    }

    static void for_each_copy_array(auto& tile_particles, auto&& fn)
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    {
        for (auto& cell : tile_particles()) // all cells, for levelghosts, batch spans cells
            fn(cell);
    }

    static void for_each_copy_array(auto& tile_particles, auto&& fn)
        requires(Particles_t::layout_mode == LayoutMode::AoSCMTS)
    {
        fn(tile_particles); // one flat array per tile
    }

    void one_copy_tile(std::size_t const tile_idx, auto& boxings, auto& view, auto& pop)
        requires(not any_in(Particles_t::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
    {
        auto const& patch_id      = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);
        auto const& patch_box     = pps.box();

        auto& tile           = pps()[tile_idx];
        bool const is_border = patch_box * grow(tile, 1) != grow(tile, 1);
        auto const& layout   = electromag.E[0][tile_idx].layout();
        auto const tile_em   = em_tile(tile_idx);
        auto& rhoP           = pop.particleDensity()[tile_idx];
        auto F               = tile_at(pop.flux(), tile_idx);
        Interpolator_t interp;

        auto const tile_cell = pps.local_cell(tile.lower);

        auto on_tile = [&]<bool border>() {
            for (auto p : tile())
            {
                per_particle(p, layout, tile_cell, 0, tile_em);
                if constexpr (!border)
                    interp.particleToMesh(p, rhoP(), F, layout);
                else
                {
                    if (isIn(p.iCell(), patch_boxings.nonLevelGhostBox))
                        interp.particleToMesh(p, rhoP(), F, layout);
                }
            }
        };

        if (is_border)
            on_tile.template operator()<true>();
        else
            on_tile.template operator()<false>();
    }

    void per_copy_of_cpu_tile(auto& boxings, auto& view, auto& pop)
    {
        for (std::size_t tile_idx = 0; tile_idx < pps().size(); ++tile_idx)
            one_copy_tile(tile_idx, boxings, view, pop);
    }

    auto em_tile(auto const tidx) _PHARE_ALL_FN_
    {
        return electromag.template as<Electromag_vt>([&] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>([&](auto i) { return vf[i][tidx](); });
        });
    }


    ParticleArray_v pps;
    Electromag_t::Super const electromag;
    double const coef1;
    std::array<double, dim> halfdt;
    Particles_t& particles;
};

} // namespace PHARE::core


#if PHARE_HAVE_MKN_GPU

namespace PHARE::core::mkn_xyz
{

// ── mkn.gpu dispatch: one ThreadedStreamLauncher host thread + stream per patch ────────
// same per patch steps/functors as core::MultiBoris, plus the GPU kernel paths

template<typename ModelAccessor, typename Interpolator>
struct MultiBoris : core::MultiBoris<ModelAccessor, Interpolator>
{
    using Super           = core::MultiBoris<ModelAccessor, Interpolator>;
    using Backend         = Super::Backend;
    using ParticleArray_t = Super::ParticleArray_t;
    using Box_t           = Box<int, Super::dim>;
    using StreamLauncher  = gpu::ThreadedStreamLauncher<ModelAccessor>;
    using GpuBoxSpanSet_t = SpanSet<Box_t, default_span_size_t, mkn::gpu::ManagedAllocator<Box_t>>;

    static constexpr auto opts = MultiBorisOptions{};

    MultiBoris(double const dt_, ModelAccessor& _accessor)
        : Super{dt_, _accessor}
    {
    }

    template<MultiBorisMode mode = MultiBorisMode::REF>
    void move(auto const& boxings);

    StreamLauncher streamer{this->accessor, opts.use_main_thread ? 0 : 1};
    GpuBoxSpanSet_t gpu_nlgb;
    GpuBoxSpanSet_t gpu_domain;
    GpuBoxSpanSet_t gpu_non_owned;
};


template<typename ModelAccessor, typename Interpolator>
template<MultiBorisMode mode>
void MultiBoris<ModelAccessor, Interpolator>::move(auto const& boxings)
{
    static constexpr auto copy     = mode == MultiBorisMode::COPY;
    static constexpr auto is_cpu   = ParticleArray_t::alloc_mode == AllocatorMode::CPU;
    static constexpr auto is_gpu   = ParticleArray_t::alloc_mode == AllocatorMode::GPU_UNIFIED;
    static constexpr auto gpu_copy = copy and is_gpu and Backend::has_gpu_copy;

    if constexpr (gpu_copy)
        Backend::prepare_gpu_copy(*this, boxings);

    auto move = [&](auto const i) mutable {
        if constexpr (gpu_copy)
            Backend::template move_gpu_copy<mode>(*this, i);
        else if constexpr (copy) // CPU, or GPU_UNIFIED without a GPU copy kernel (AoSPCTS)
        {
            if constexpr (is_gpu) // unified memory: wait for prior work on this patch
                streamer.streams[i].sync();
            Backend::template move_cpu_copy<mode>(*this, boxings, i);
        }
        else
            Backend::template move_rest<mode>(*this, i);
    };
    auto sync = [&](auto const i) mutable {
        if constexpr (not copy)
            Backend::sync_ref(*this, i);
    };

    if constexpr (opts.use_main_thread)
    {
        for (std::size_t i = 0; i < this->accessor.size(); ++i)
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
        streamer.host(std::move(move));
        streamer.host(std::move(sync));
    }
}


} // namespace PHARE::core::mkn_xyz

#endif // PHARE_HAVE_MKN_GPU


#endif
