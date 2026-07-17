#ifndef PHARE_BENCH_CORE_INTERPOLATOR
#define PHARE_BENCH_CORE_INTERPOLATOR

#include "phare_core.hpp"
#include "bench/core/bench.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

template<std::size_t dim, std::size_t interp, typename ParticleArray_t, bool sort_particles = false>
void interpolate(benchmark::State& state)
{
    using namespace PHARE::core;

    std::uint32_t constexpr static cells   = 30;
    std::uint32_t constexpr static n_parts = 1e6;
    auto constexpr static alloc_mode       = ParticleArray_t::alloc_mode;
    bool constexpr static c_ordering       = true;

    using PHARE_Types     = PHARE_Types<dim, interp>;
    using GridLayout_t    = typename PHARE_Types::GridLayout_t;
    using Field_t         = typename PHARE_Types::Field_t;
    using NdArrayVector_t = NdArrayVector<dim, double, c_ordering, alloc_mode>;
    using Grid_t          = Grid<NdArrayVector_t, HybridQuantity::Scalar>;
    using Electromag_t    = UsableElectromag<GridLayout_t, ParticleArray_t::alloc_mode>;
    using VecField_t      = UsableVecField<dim, ParticleArray_t::alloc_mode>;

    TestGridLayout<GridLayout_t> layout{cells};
    Electromag_t em{layout};
    VecField_t flux{"F", layout, HybridQuantity::Vector::V};

    ParticleArray_t particles{n_parts, particle<dim>()};
    disperse(particles, layout.AMRBox(), 133337);

    Grid_t rho{"rho", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho)};

    if constexpr (sort_particles)
        ParticleArraySorter<ParticleArray_t>{particles, layout.AMRBox()}();

    if constexpr (alloc_mode == PHARE::AllocatorMode::CPU)
    {
        Interpolator<dim, interp> interpolator;
        while (state.KeepRunningBatch(1))
        {
            interpolator(particles, rho, flux, layout);
        }
    }
    else if constexpr (alloc_mode == PHARE::AllocatorMode::GPU_UNIFIED)
    {
        // take views for GPU lambda copy
        auto em_view   = *em;
        auto flux_view = *flux;
        auto rho_view  = *rho;
        auto ps        = *particles;

        while (state.KeepRunningBatch(1))
        {
            constexpr static bool atomic_ops = true;
            PHARE_WITH_MKN_GPU( //
                mkn::gpu::GDLauncher{n_parts}([=] __device__() mutable {
                    Interpolator<dim, interp, atomic_ops> interpolator;
                    auto particle = ps.begin() + mkn::gpu::idx();
                    interpolator.particleToMesh(deref(particle), rho_view, flux_view, layout);
                }); //
            )
        }
    }
    else
        std::abort();
}

#endif /*PHARE_BENCH_CORE_INTERPOLATOR*/
