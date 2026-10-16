#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP

#include "core/utilities/kernels.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/utilities/range/range.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "interpolator.hpp"


namespace PHARE::core
{

// visits, for tile tidx, the particles of every cell its level ghost halo reaches, read
// from that cell's clamp-owner tile (TileSet::tag_cells_). Level ghosts are never in the
// patch domain, so only halo cells outside of it are visited: none for interior tiles
template<typename Particles>
void on_reachable_level_ghosts(Particles const& particles, std::size_t const tidx, auto&& fn)
    requires(any_in(Particles::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
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


// simple facade to launch e.g. GPU kernels if needed
//  and to not need to modify the interpolator much for gpu specifically

template<std::size_t dim, std::size_t interpOrder, bool atomic_ops = false,
         typename Interpolator_t = Interpolator<dim, interpOrder, atomic_ops>>
class Interpolating
{
public:
    template<typename Particles>
    inline void operator()(Particles const& particles, auto& rhoP, auto& rhoC, auto& flux,
                           auto const& layout, double coef = 1.)
    {
        PHARE_LOG_SCOPE(3, "Interpolating::particleToMesh");
        particleToMesh(particles, layout, rhoP, rhoC, flux, coef);
    }

    // Tiled LevelGhost: level ghosts live only in their clamp-owner tile, and their deposits
    // happen after the tile reduction, so each tile gathers (read-only) every level ghost
    // its halo reaches from the owning tile.
    // GPU_UNIFIED falls back to the CPU path until there's an optimized GPU version

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(any_in(Particles::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
    _PHARE_ALL_FN_
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
            particleToMesh_tile<type>(tidx, particles, rhoP, rhoC, flux, coef);
    }

    // only writes the moments of tile tidx, distinct tiles can be deposited concurrently
    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh_tile(std::size_t const tidx, Particles const& particles, auto& rhoP,
                             auto& rhoC, auto& flux, double coef = 1.)
        requires(any_in(Particles::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS))
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        auto [rhop, rhoc, F] = tiles_at(tidx, rhoP, rhoC, flux);
        auto const& rho_lay  = tile_layout(rhoP, tidx);
        auto const deposit   = [&](auto const& cell_particles) {
            for (auto const& p : cell_particles)
                interp_.particleToMesh(p, rhop, rhoc, F, rho_lay, coef);
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
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSMapped)
    _PHARE_ALL_FN_
    {
        interp_(particles, rhoP, rhoC, flux, layout, coef);
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSTS)
    _PHARE_ALL_FN_
    {
        // GPU_UNIFIED falls back to the CPU path until there's an optimized GPU version
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto [rhop, rhoc, F] = tiles_at(tidx, rhoP, rhoC, flux);
            auto const& rho_lay  = tile_layout(rhoP, tidx);
            for (auto const& p : particles()[tidx]())
                interp_.particleToMesh(p, rhop, rhoc, F, rho_lay, coef);
        }
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::SoATS)
    _PHARE_ALL_FN_
    {
        for (auto const& tile : particles())
            for (std::size_t i = 0; i < tile().size(); ++i)
                interp_.particleToMesh(tile()[i], rhoP, rhoC, flux, layout, coef);
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(Particles::layout_mode == LayoutMode::AoSPC)
    _PHARE_ALL_FN_
    {
        auto constexpr static alloc_mode = Particles::alloc_mode;

        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            for (auto const& bix : particles.local_box())
                for (auto const& p : particles(bix))
                    interp_.particleToMesh(p, rhoP, rhoC, flux, layout, coef);
        }
        // else if (alloc_mode == AllocatorMode::GPU_UNIFIED and impl < 2)
        // {
        //     static_assert(atomic_ops, "GPU must be atomic");
        //     PHARE_WITH_MKN_GPU(
        //         mkn::gpu::GDLauncher{particles.size()}([=] _PHARE_ALL_FN_() mutable {
        //             auto it = particles[mkn::gpu::idx()];
        //             Interpolator_t{}.particleToMesh(*it, density, flux, layout, coef);
        //         }); //
        //     )
        // }
        else if (alloc_mode == AllocatorMode::GPU_UNIFIED /*and impl == 2*/)
        {
            PHARE_WITH_MKN_GPU({ //
                using lobox_t  = Particles::lobox_t;
                using Launcher = gpu::BoxCellNLauncher<lobox_t>;
                auto lobox     = particles.local_box();
                Launcher{lobox, particles.max_size()}([=] _PHARE_ALL_FN_() mutable {
                    box_kernel(particles, layout, flux, rhoP, rhoC);
                });
            })
        }
        else
            throw std::runtime_error("fail");
    }

    template<auto type = ParticleType::Domain, typename Particles>
    void particleToMesh(Particles const& particles, auto const& layout, auto& rhoP, auto& rhoC,
                        auto& flux, double coef = 1.)
        requires(any_in(Particles::layout_mode, LayoutMode::AoS, LayoutMode::SoA,
                        LayoutMode::SoAVX))
    _PHARE_ALL_FN_
    {
        auto constexpr static alloc_mode = Particles::alloc_mode;

        if constexpr (alloc_mode == AllocatorMode::CPU)
        {
            if constexpr (Particles::layout_mode == LayoutMode::SoAVX)
            {
                for (auto const& it : particles)
                    interp_.particleToMesh(it, rhoP, rhoC, flux, layout, coef);
            }
            else
                interp_(particles, rhoP, rhoC, flux, layout, coef);
        }
        else if (alloc_mode == AllocatorMode::GPU_UNIFIED)
        {
            static_assert(atomic_ops, "GPU must be atomic");
            PHARE_WITH_MKN_GPU( //
                kernel::launch(particles.size(),
                               [=] _PHARE_ALL_FN_() mutable {
                                   auto particle = particles.begin() + kernel::idx();
                                   Interpolator_t{}.particleToMesh(deref(particle), rhoP, rhoC,
                                                                   flux, layout, coef);
                               }); //
            )
        }
        else
            throw std::runtime_error("fail");
    }

    // AoSPCTS GPU deposit, launched with gpu::TileCellLauncher: block = (tile, cell), the block's
    // threads stride the cell's particles. Cells share field nodes, so the interpolator must be
    // atomic. Field_t as on_tiles (tiled field value_type)
    template<typename Field_t, typename Launcher>
    static void on_tile_cells(auto& particles, auto& Flux, auto& Rhop, auto& Rhoc,
                              std::size_t const threads, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        static_assert(atomic_ops, "GPU must be atomic");
        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;

        auto const tidx   = Launcher::tile_idx();
        auto& cps         = particles()[tidx]();
        auto const& lobox = cps.local_box();
        if (Launcher::cell_idx() >= lobox.size())
            return;

        auto& parts         = cps(*(lobox.begin() + Launcher::cell_idx()));
        auto& rhop_tile     = Rhop[tidx];
        auto& rhoc_tile     = Rhoc[tidx];
        auto flux           = Flux.template as<TensorField_t>([&](auto& c) { return c[tidx](); });
        auto const& layout  = rhop_tile.layout();
        auto& pdensity      = rhop_tile();
        auto& cdensity      = rhoc_tile();
        for (std::size_t pid = Launcher::thread_idx(); pid < parts.size(); pid += threads)
            Interpolator_t{}.particleToMesh(parts[pid], pdensity, cdensity, flux, layout, coef);
#endif // PHARE_HAVE_MKN_GPU_HW
    }

    template<bool in_box = false, typename Particles, typename GridLayout, typename VecField,
             typename Field>
    static void box_kernel(Particles& particles, GridLayout& layout, VecField& flux, Field& density,
                           double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        using lobox_t         = Particles::lobox_t;
        using Launcher        = gpu::BoxCellNLauncher<lobox_t>;
        auto const& lobox     = particles.local_box();
        auto const& blockidx  = Launcher::block_idx();
        auto const& threadIdx = Launcher::thread_idx();
        auto& parts           = particles(*(lobox.begin() + blockidx));

        auto const interp = [&]() {
            Interpolator_t{}.particleToMesh(parts[threadIdx], density, flux, layout, coef);
        };

        if (threadIdx < parts.size())
        {
            if constexpr (!in_box)
            {
                interp();
            }
            else
            {
                if (isIn(parts[threadIdx], particles.box()))
                    interp();
            }
        }
#endif // PHARE_HAVE_MKN_GPU_HW
    }

#if 0 // old unused
    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel_ts(Particles& particles, GridLayout& layout, VecField& flux,
                                Field& density, double coef = 1.) _PHARE_ALL_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW

        Point<std::uint32_t, 3> const rcell{threadIdx.z, threadIdx.y, threadIdx.x};
        extern __shared__ double data[];

        auto& tile        = *particles().at(rcell + 2);
        auto& parts       = tile();
        auto const each   = parts.size() / 8;
        auto const tile_p = rcell % 2;
        auto const pidx   = tile_p[2] + tile_p[1] * 2 + tile_p[0] * 2 * 2;
        auto const ji     = [&](auto const i) { return i * 8 + pidx; };
        auto doX          = [&]<std::uint8_t IDX = 0>(auto& feeld) {
            auto v = make_array_view(&data[0], feeld.shape());
            if (kernel::idx() == 0)
                for (auto& e : v)
                    e = 0;
            Interpolator_t interp;
            std::size_t i = 0, j = ji(i);
            __syncthreads();
            for (; i < each; ++i, j = ji(i))
                interp.p2m_setup(parts[j], layout),
                    interp.template p2m_per_component<IDX>(parts[j], v);
            if (pidx < parts.size() - (8 * each))
                interp.p2m_setup(parts[j], layout),
                    interp.template p2m_per_component<IDX>(parts[j], v);
            __syncthreads();
            if (kernel::idx() == 0)
                for (i = 0; i < v.size(); ++i)
                    feeld.data()[i] = v.data()[i];
        };

        doX(density);
        doX.template operator()<1>(flux[0]);
        doX.template operator()<2>(flux[1]);
        doX.template operator()<3>(flux[2]);


#endif // PHARE_HAVE_MKN_GPU_HW
    }
#endif

    template<typename Field_t>
    static void on_tiles(auto& particles, auto& Flux, auto& Rhop, auto& Rhoc,
                         double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        // printf("L:%d i %llu \n", __LINE__, blockIdx.x);
        static_assert(PHARE_HAVE_MKN_GPU_HW);

        static_assert(atomic_ops, "GPU must be atomic");
        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;

        extern __shared__ double data[];

        auto& rhop_tile = Rhop[blockIdx.x];
        auto& rhoc_tile = Rhoc[blockIdx.x];
        auto ft         = Flux.template as<TensorField_t>([&](auto& c) { return c[blockIdx.x](); });
        auto& pdensity  = rhop_tile(); // as ref = crash?
        auto& cdensity  = rhoc_tile(); // as ref = crash?

        auto& tile                        = particles()[blockIdx.x];
        auto const tbox                   = rhop_tile.ghost_box();
        auto const ziz                    = product(tbox.shape());
        auto const tidx                   = threadIdx.x;
        auto const ws                     = kernel::warp_size();
        auto const& [flux0, flux1, flux2] = ft();
        auto const r0                     = pdensity;
        auto const r1                     = cdensity;
        auto flux                         = ft;

        auto rhop = make_array_view(&data[ziz * 0], pdensity.shape());
        auto rhoc = make_array_view(&data[ziz * 1], cdensity.shape());
        auto fx   = make_array_view(&data[ziz * 2], flux[0].shape());
        auto fy   = make_array_view(&data[ziz * 3], flux[1].shape());
        auto fz   = make_array_view(&data[ziz * 4], flux[2].shape());

        pdensity.reset(rhop);
        cdensity.reset(rhoc);
        flux[0].reset(fx);
        flux[1].reset(fy);
        flux[2].reset(fz);

        for (int i = tidx; i < ziz * 5; i += blockDim.x)
            data[i] = 0;


        auto& parts = tile();
        {
            Interpolator_t interp;

            __syncthreads();

            for (std::size_t i = tidx; i < parts.size(); i += blockDim.x)
                interp.particleToMesh(parts[i], pdensity, cdensity, flux, tbox, coef);

            __syncthreads();
        }


        pdensity.reset(r0);
        cdensity.reset(r1);

        for (int i = tidx; i < ziz; i += blockDim.x)
        {
            pdensity.data()[i] += rhop.data()[i];
            cdensity.data()[i] += rhoc.data()[i];
            flux0.data()[i] += fx.data()[i];
            flux1.data()[i] += fy.data()[i];
            flux2.data()[i] += fz.data()[i];
        }

#endif
    }


    // Combined Boris push + particle-to-mesh deposit in a single GPU kernel.
    // Originals are untouched: each thread works on a register copy.
    // Qualifying particles (isIn filter_box and not in exclude_boxes after push) are atomically
    // placed into a shared warp buffer. Once the buffer holds a full warp the whole warp
    // deposits together; any remainder is flushed at the end.
    //
    // Shared memory layout (launcher.ds must cover):
    //   [5 * ziz * sizeof(double)]            field accumulators
    //   [2 * warp_size * sizeof(Particle_t)]  warp buffer
    //   [sizeof(int)]                          atomic fill counter
    //
    // 2*warp_size slots are the minimum needed: before a batch buf_count can be
    // warp_size-1, and if every thread in the batch qualifies buf_count reaches
    // 2*warp_size-1 — still within bounds.
    template<typename Field_t, auto particle_type_v, auto alloc_mode_v>
    static void on_tiles_push_deposit(auto& particles, auto const& electromag, auto& Flux,
                                      auto& Rhop, auto& Rhoc, auto const& filter_box,
                                      auto const& exclude_boxes, double const dto2m,
                                      auto const& halfdt,
                                      double const coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        static_assert(atomic_ops, "GPU must be atomic");

        using TensorField_t = basic::TensorField<typename Field_t::value_type, 1>;
        using Scalar_vt     = typename Field_t::value_type;
        using VecField_vt   = basic::TensorField<Scalar_vt, 1>;
        using Electromag_vt = basic::Electromag<VecField_vt>;
        using Particle_t    = ParticleDefaults<dim>::Particle_t;

        extern __shared__ double data[];

        auto& rhop_tile = Rhop[blockIdx.x];
        auto& rhoc_tile = Rhoc[blockIdx.x];
        auto ft         = Flux.template as<TensorField_t>([&](auto& c) { return c[blockIdx.x](); });
        auto& pdensity  = rhop_tile();
        auto& cdensity  = rhoc_tile();

        auto const tbox                   = rhop_tile.ghost_box();
        auto const ziz                    = product(tbox.shape());
        auto const tidx                   = threadIdx.x;
        auto const ws                     = kernel::warp_size();
        auto const& [flux0, flux1, flux2] = ft();
        auto const r0                     = pdensity;
        auto const r1                     = cdensity;
        auto flux                         = ft;

        auto rhop = make_array_view(&data[ziz * 0], pdensity.shape());
        auto rhoc = make_array_view(&data[ziz * 1], cdensity.shape());
        auto fx   = make_array_view(&data[ziz * 2], flux[0].shape());
        auto fy   = make_array_view(&data[ziz * 3], flux[1].shape());
        auto fz   = make_array_view(&data[ziz * 4], flux[2].shape());

        pdensity.reset(rhop);
        cdensity.reset(rhoc);
        flux[0].reset(fx);
        flux[1].reset(fy);
        flux[2].reset(fz);

        // Warp buffer and atomic counter after the field region.
        auto* const part_buf = reinterpret_cast<Particle_t*>(reinterpret_cast<char*>(&data[0])
                                                             + ziz * 5 * sizeof(double));
        auto& buf_count      = *reinterpret_cast<int*>(reinterpret_cast<char*>(part_buf)
                                                       + 2 * ws * sizeof(Particle_t));

        for (int i = tidx; i < static_cast<int>(ziz * 5); i += blockDim.x)
            data[i] = 0;
        if (tidx == 0)
            buf_count = 0;

        __syncthreads();

        // Per-tile EM view and layout
        auto const em = electromag.template as<Electromag_vt>([=] _PHARE_ALL_FN_(auto const& vf) {
            return for_N_make_array<3>(
                [&vf] _PHARE_ALL_FN_(auto i) { return vf[i][blockIdx.x](); });
        });
        auto const& layout = electromag.E[0][blockIdx.x].layout();

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
                    boris::accelerate(p, interp.m2p(p, em, layout), dto2m);
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
                interp.particleToMesh(part_buf[tidx], pdensity, cdensity, flux, tbox, coef);

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
            interp.particleToMesh(part_buf[tidx], pdensity, cdensity, flux, tbox, coef);
        }

        __syncthreads();

        pdensity.reset(r0);
        cdensity.reset(r1);

        for (int i = tidx; i < static_cast<int>(ziz); i += blockDim.x)
        {
            pdensity.data()[i] += rhop.data()[i];
            cdensity.data()[i] += rhoc.data()[i];
            flux0.data()[i] += fx.data()[i];
            flux1.data()[i] += fy.data()[i];
            flux2.data()[i] += fz.data()[i];
        }
#endif
    }


    template<typename Field, typename Particles, typename GridLayout>
    static void on_cpu_tiles(Particles& particles, GridLayout& layout,
                             double coef = 1.) _PHARE_ALL_FN_
    {
        using TensorField_t = basic::TensorField<typename Field::Super, 1>;

        for (auto& tile : particles())
        {
            auto& parts                                = tile();
            auto const tbox                            = tile.field_box();
            auto const& [density, flux0, flux1, flux2] = tile.fields();
            auto const r0                              = density;

            TensorField_t flux{flux0, flux1, flux2};
            Interpolator_t interp;
            std::size_t pid = 0;

            interp.particleToMesh(parts[pid], density, flux, /*layout,*/ tbox, coef);
        }
    }




    template<typename Particles, typename GridLayout, typename VecField, typename Field>
    static void chunk_kernel(Particles& particles, GridLayout& layout, VecField& flux,
                             Field& density, double coef = 1.) _PHARE_DEV_FN_
    {
#if PHARE_HAVE_MKN_GPU_HW
        auto constexpr static nghosts = GridLayout::options.field_ghost_width;
        extern __shared__ double data[];
        auto const& lobox = particles.local_box();
        auto const ziz    = 9 * 9 * 9;

        auto const t_x = threadIdx.x;
        auto const t_y = threadIdx.y;
        auto const t_z = threadIdx.z;

        core::Point<std::uint32_t, 3> const tcell{t_x, t_y, t_z};
        core::Point<std::uint32_t, 3> const locell = tcell + nghosts;

        {
            Interpolator_t interp;
            auto& parts = particles(locell);

            { // rho
                auto const r0 = density.data();
                auto v        = make_array_view(&data[0], density.shape());
                density.setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], v);
                }
                __syncthreads();
                density.setData(r0);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        density.data()[i] = v.data()[i];
            }

            {
                auto const f0 = flux[0].data();
                auto v        = make_array_view(&data[0], flux[0].shape());
                flux[0].setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[0]);
                }
                __syncthreads();
                flux[0].setData(f0);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[0].data()[i] = v.data()[i];
            }

            {
                auto const f1 = flux[1].data();
                auto v        = make_array_view(&data[0], flux[1].shape());
                flux[1].setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[1]);
                }
                __syncthreads();
                flux[1].setData(f1);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[1].data()[i] = v.data()[i];
            }

            {
                auto const f2 = flux[2].data();
                auto v        = make_array_view(&data[0], flux[2].shape());
                flux[2].setData(v.data());
                if (kernel::idx() == 0)
                    for (auto& e : v)
                        e = 0;
                __syncthreads();
                for (std::size_t i = 0; i < parts.size(); ++i)
                {
                    p2m_setup(parts[i], layout);
                    p2m_per_component(parts[i], flux[2]);
                }
                __syncthreads();
                flux[2].setData(f2);
                if (kernel::idx() == 0)
                    for (std::size_t i = 0; i < ziz; ++i)
                        flux[2].data()[i] = v.data()[i];
            }
        }


#endif // PHARE_HAVE_MKN_GPU_HW
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
        assert(momentumTensor.isUsable());
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
        requires(Particles_t::layout_mode == LayoutMode::AoSPCTS)
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt            = tile_at(momentumTensor, tidx);
            auto const& mt_lay = tile_layout(momentumTensor, tidx);
            auto const deposit
                = [&](auto const& cell_particles) { interp_(cell_particles, mt, mt_lay, mass); };

            if constexpr (type == ParticleType::LevelGhost)
                on_reachable_level_ghosts(particles, tidx, deposit);
            else
            {
                auto& cps = particles()[tidx]();
                for (auto const& bix : cps.local_box())
                    deposit(cps(bix));
            }
        }
    }

    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSCMTS)
    {
        static_assert(any_in(type, ParticleType::Domain, ParticleType::LevelGhost));

        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
        {
            auto mt            = tile_at(momentumTensor, tidx);
            auto const& mt_lay = tile_layout(momentumTensor, tidx);
            auto const deposit
                = [&](auto const& tile_particles) { interp_(tile_particles, mt, mt_lay, mass); };

            if constexpr (type == ParticleType::LevelGhost)
                on_reachable_level_ghosts(particles, tidx, deposit);
            else
                deposit(particles()[tidx]());
        }
    }

    template<auto type = ParticleType::Domain, typename Particles_t>
    void particleToMesh(Particles_t& particles, auto& momentumTensor, auto const& layout,
                        double mass = 1.)
        requires(Particles_t::layout_mode == LayoutMode::AoSTS)
    {
        for (std::size_t tidx = 0; tidx < particles().size(); ++tidx)
            interp_(particles()[tidx](), tile_at(momentumTensor, tidx),
                    tile_layout(momentumTensor, tidx), mass);
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATING_HPP*/
