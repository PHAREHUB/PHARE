// IWYU pragma: private, include "core/numerics/pusher/any_boris.hpp"

#ifndef PHARE_CORE_PUSHER_BORIS_MKN_ANY_BORIS_HPP
#define PHARE_CORE_PUSHER_BORIS_MKN_ANY_BORIS_HPP

#include "core/utilities/kernels.hpp"
#include "core/numerics/pusher/boris/basics.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_service.hpp"


namespace PHARE::core
{

// Primary template — generic CPU/GPU fallback for non-specialized layouts.
// Specializations in mkn_any_boris_aos.hpp override for specific (layout, alloc) pairs.
template<LayoutMode layout, AllocatorMode alloc, std::size_t dim, typename Interpolator,
         typename GridLayout>
struct AnyBorisImpl
{
    template<auto alloc_m, typename P, typename E, typename I, typename G>
    static void pp(P& particle, E const& em, I& interp, G const& gl,
                   std::array<double, dim> const& halfDtOverDl, double const dto2m) _PHARE_ALL_FN_
    {
        particle.iCell() = boris::advance<alloc_m>(particle, halfDtOverDl);
        boris::accelerate(particle, interp.m2p(particle, em, gl), dto2m);
        particle.iCell() = boris::advance<alloc_m>(particle, halfDtOverDl);
    }

    template<auto type, typename Particles, typename Electromag>
    static void move(Particles& particles, Electromag const& em, Interpolator& interp,
                     GridLayout const& gl, std::array<double, dim> const& halfDtOverDl,
                     double const dto2m) _PHARE_ALL_FN_
    {
        if constexpr (alloc == AllocatorMode::CPU)
        {
            if constexpr (layout == LayoutMode::SoAVX)
                boris::avx_ad_acc_ad(particles, em, interp, gl, halfDtOverDl, dto2m);
            else
                for (auto& p : particles)
                    pp<alloc>(p, em, interp, gl, halfDtOverDl, dto2m);
        }
        else if constexpr (alloc == AllocatorMode::GPU_UNIFIED)
        {
            auto view = *particles;
            PHARE_ASSERT(view.size());
            kernel::launch(particles.size(), [=, gl_ = gl, halfDtOverDl_ = halfDtOverDl,
                                              dto2m_ = dto2m] _PHARE_ALL_FN_() mutable {
                Interpolator interp_d;
                auto particle = view.begin() + kernel::idx();
                pp<alloc>(particle, em, interp_d, gl_, halfDtOverDl_, dto2m_);
            });
        }
        else
            throw std::runtime_error("AnyBorisImpl: unsupported alloc mode");
    }
};


template<std::size_t dim, typename Interpolator, typename GridLayout>
struct AnyBorisImpl<LayoutMode::AoSPC, AllocatorMode::CPU, dim, Interpolator, GridLayout>
{
    template<auto type, typename Particles, typename Electromag>
    static void move(Particles& particles, Electromag const& em, Interpolator& interp,
                     GridLayout const& gl, std::array<double, dim> const& halfDtOverDl,
                     double const dto2m)
    {
        PHARE_LOG_SCOPE(1, "boris::move_aos_pc_cpu");
        constexpr auto alloc_mode = AllocatorMode::CPU;

        auto per_particle = [&]<bool accelerate = false>() mutable {
            for (auto const& bix : particles.local_box())
                for (std::size_t idx = 0; idx < particles.size(bix); ++idx)
                {
                    auto& particle = particles(bix)[idx];
                    if constexpr (accelerate)
                        boris::accelerate(particle, interp.m2p(particle, em, gl), dto2m);
                    auto const& newCell = boris::advance<alloc_mode>(particle, halfDtOverDl);
                    if (newCell != particle.iCell())
                    {
                        ParticleTracker<dim> const pt{particle.iCell()};
                        particle.iCell() = newCell;
                        particles.template move_check<type>(pt, idx, particle);
                    }
                }
        };

        {
            PHARE_LOG_SCOPE(1, "boris::move_aos_pc_cpu::advance_0");
            per_particle();
        }
        ParticleArrayService::sync<0, type>(particles);
        {
            PHARE_LOG_SCOPE(1, "boris::move_aos_pc_cpu::advance_accelerate");
            per_particle.template operator()<true>();
        }
        ParticleArrayService::sync<1, type>(particles);
    }
};


template<std::size_t dim, typename Interpolator, typename GridLayout>
struct AnyBorisImpl<LayoutMode::AoSPC, AllocatorMode::GPU_UNIFIED, dim, Interpolator, GridLayout>
{
    template<auto type, typename Particles, typename Electromag>
    static void move(Particles& particles, Electromag const& em, Interpolator& /*interp*/,
                     GridLayout const& gl, std::array<double, dim> const& halfDtOverDl,
                     double const dto2m)
    {
        PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu");

        if constexpr (any_in(Particles::impl_v, 0, 1))
            move_impl<type>(particles, em, gl, halfDtOverDl, dto2m);
        else if constexpr (any_in(Particles::impl_v, 2))
            move_impl_2<type>(particles, em, gl, halfDtOverDl, dto2m);
        else
            throw std::runtime_error("Unhandled AoSPC Boris GPU impl");
    }

    template<auto type, typename Particles, typename Electromag>
    static void move_impl(Particles& particles, Electromag const& em, GridLayout const& gl,
                          std::array<double, dim> const& halfDtOverDl, double const dto2m)
    {
        PHARE_WITH_MKN_GPU({
            constexpr auto alloc_mode = AllocatorMode::GPU_UNIFIED;
            auto view                 = *particles;

            auto per_particle =
                [=, gl_ = gl, halfDtOverDl_ = halfDtOverDl, dto2m_ =
                 dto2m] _PHARE_ALL_FN_<bool accelerate = false>() mutable
            {
                Interpolator interp;
                auto it        = view[kernel::idx()];
                auto& particle = *it;
                if constexpr (accelerate)
                    boris::accelerate(particle, interp.m2p(particle, em, gl_), dto2m_);
                auto const& newCell = boris::advance<alloc_mode>(particle, halfDtOverDl_);
                if (!array_equals(newCell, particle.iCell()))
                {
                    it.icell_changer(newCell);
                    particle.iCell() = newCell;
                }
            };

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_0/1");
                kernel::launch(particles.size(), per_particle);
            }
            if constexpr (any_in(Particles::impl_v, 1))
                ParticleArrayService::sync<0, type>(view);
            ParticleArrayService::sync<0, type>(particles);

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_accelerate_0/1");
                kernel::launch(particles.size(), [=] _PHARE_ALL_FN_() mutable {
                    per_particle.template operator()<true>();
                });
            }
            if constexpr (any_in(Particles::impl_v, 1))
                ParticleArrayService::sync<1, type>(view);
            ParticleArrayService::sync<1, type>(particles);
        })
    }

    template<auto type, typename Particles, typename Electromag>
    static void move_impl_2(Particles& particles, Electromag const& em, GridLayout const& gl,
                            std::array<double, dim> const& halfDtOverDl, double const dto2m)
    {
        PHARE_WITH_GPU({
            constexpr auto alloc_mode = AllocatorMode::GPU_UNIFIED;
            auto view                 = *particles;

            auto per_particle =
                [=, gl_ = gl, halfDtOverDl_ = halfDtOverDl, dto2m_ =
                 dto2m] _PHARE_ALL_FN_<bool accelerate = false>() mutable
            {
                auto const blockidx  = kernel::block_idx();
                auto const threadIdx = kernel::thread_idx();
                auto const& locell   = *(view.local_box().begin() + blockidx);
                auto& parts          = view(locell);
                if (threadIdx >= parts.size())
                    return;
                auto& particle = parts[threadIdx];
                Interpolator interp;
                if constexpr (accelerate)
                    boris::accelerate(particle, interp.m2p(particle, em, gl_), dto2m_);
                auto const& newCell = boris::advance<alloc_mode>(particle, halfDtOverDl_);
                if (!array_equals(newCell, particle.iCell()))
                {
                    ParticleTracker<dim> const pt{particle.iCell()};
                    particle.iCell() = newCell;
                    view.template move_check<type>(pt, threadIdx, particle);
                }
            };

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_2");
                kernel::launch(view.local_box(), particles.max_size(), per_particle);
            }
            ParticleArrayService::sync<0, type>(view);
            ParticleArrayService::sync<0, type>(particles);

            {
                PHARE_LOG_SCOPE(1, "boris::move_aos_pc_gpu::advance_accelerate_2");
                kernel::launch(
                    view.local_box(), particles.max_size(),
                    [=] _PHARE_ALL_FN_() mutable { per_particle.template operator()<true>(); });
            }
            ParticleArrayService::sync<1, type>(view);
            ParticleArrayService::sync<1, type>(particles);
        })
    }
};

} // namespace PHARE::core


#endif
