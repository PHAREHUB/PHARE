#ifndef PHARE_SOLVER_GPU_OFFLOADER_H
#define PHARE_SOLVER_GPU_OFFLOADER_H

#if !defined(PHARE_WITH_GPU)
#error "PHARE_WITH_GPU undefined"
#endif

#include "core/utilities/types.h"

#define PHARE_WITH_GPU_MKN // TORM

namespace PHARE::solver::gpu
{
template<typename Solver>
class BaseOffloader
{
public:
    static constexpr auto dimension    = Solver::dimension;
    static constexpr auto interp_order = Solver::interp_order;

    using GridLayout = typename Solver::GridLayout;

    virtual ~BaseOffloader() = default;
};
} // namespace PHARE::solver::gpu

#if defined(PHARE_WITH_GPU_MKN)
#include "kul/gpu/asio.hpp"
#include "core/numerics/pusher/granov.h"
namespace PHARE::solver::gpu_mkn
{
template<typename GridLayout_>
struct PatchState
{
    using Float                     = double;
    using GridLayout                = GridLayout_;
    static constexpr auto dimension = GridLayout::dimension;
    using HybridQuantity            = core::HybridQuantity;

    template<typename State>
    PatchState(GridLayout const& gridLayout, State& state)
        : layout{gridLayout}
    {
        auto& E = state.electromag.E;
        electromag.emplace_back(view(E[0].data(), HybridQuantity::Scalar::Ex));
        electromag.emplace_back(view(E[1].data(), HybridQuantity::Scalar::Ey));
        electromag.emplace_back(view(E[2].data(), HybridQuantity::Scalar::Ez));

        auto& B = state.electromag.B;
        electromag.emplace_back(view(B[0].data(), HybridQuantity::Scalar::Bx));
        electromag.emplace_back(view(B[1].data(), HybridQuantity::Scalar::By));
        electromag.emplace_back(view(B[2].data(), HybridQuantity::Scalar::Bz));

        for (auto& pop : state.ions)
        {
            ions.emplace_back(&pop.domainParticles());
            assert(ions.back()->size()); // ?
            masses.emplace_back(pop.mass());

            density.emplace_back(view(pop.density().data(), HybridQuantity::Scalar::rho));
            assert(density.back().data() == pop.density().data());

            auto& F = pop.flux();
            flux.emplace_back(view(F[0].data(), HybridQuantity::Scalar::rho));
            flux.emplace_back(view(F[1].data(), HybridQuantity::Scalar::rho));
            flux.emplace_back(view(F[2].data(), HybridQuantity::Scalar::rho));
        }
    }

    template<typename T, typename Qty>
    static auto view(GridLayout_ const& layout, T* ptr, Qty qty)
    {
        return core::NdArrayView<dimension, T, T*>{ptr, layout.allocSize(qty)};
    }
    template<typename T, typename Qty>
    auto view(T* ptr, Qty qty)
    {
        return view(layout, ptr, qty);
    }


    GridLayout const layout;
    std::vector<core::NdArrayView<dimension, Float, Float*>> electromag;
    std::vector<core::NdArrayView<dimension, Float, Float*>> density, flux;
    std::vector<core::ParticleArray<dimension>*> ions;
    std::vector<double> masses;
};


template<std::size_t dim>
struct FieldInterop : core::NdArrayView<dim, double, double*>
{
    using Super = core::NdArrayView<dim, double, double*>;

    template<typename A, typename B>
    FieldInterop(A ptr, B nCells) __device__ : Super{ptr, nCells}
    {
    }
};

template<std::size_t dim>
struct VecFieldInterop
{
    auto __device__ getComponents() { return std::forward_as_tuple(x, y, z); }
    auto __device__ getComponents() const { return std::forward_as_tuple(x, y, z); }

    FieldInterop<dim> x, y, z;
};

template<typename GridLayout_, bool GPU = false>
struct ETC : kul::gpu::DeviceClass<GPU>
{
    using Float = double;
    using Super = kul::gpu::DeviceClass<GPU>;
    using gpu_t = ETC<GridLayout_, true>;

    template<typename T>
    using container_t = typename Super::template container_t<T>;

    ETC(PatchState<GridLayout_>& state, std::size_t pop_idx)
        : Ex{state.electromag[0].size()}
        , Ey{state.electromag[1].size()}
        , Ez{state.electromag[2].size()}
        , Bx{state.electromag[3].size()}
        , By{state.electromag[4].size()}
        , Bz{state.electromag[5].size()}
        , density{state.density[pop_idx]}
        , fluxX{state.flux[0 + (3 * pop_idx)]}
        , fluxY{state.flux[1 + (3 * pop_idx)]}
        , fluxZ{state.flux[2 + (3 * pop_idx)]}
    {
    }
    ETC(ETC const&) = delete;
    auto& operator=(ETC const&) = delete;

    auto __host__ operator()()
    {
        return Super::template alloc<gpu_t>(Ex, Ey, Ez, Bx, By, Bz, density, fluxX, fluxY, fluxZ);
    }

    struct EM_Interop
    {
        VecFieldInterop<GridLayout_::dimension> E, B;
    };


    auto __device__ operator()(GridLayout_ const& layout) const
    {
        auto constexpr dim   = GridLayout_::dimension;
        using FI             = FieldInterop<dim>;
        using VFI            = VecFieldInterop<dim>;
        using HybridQuantity = core::HybridQuantity;
        using Scalar         = HybridQuantity::Scalar;

        auto view = [&](auto ptr, auto qty) { return FI{ptr, layout.allocSize(qty)}; };

        return std::make_tuple(
            view(density, Scalar::rho),
            VFI{view(fluxX, Scalar::rho), view(fluxY, Scalar::rho), view(fluxZ, Scalar::rho)},
            EM_Interop{view(Ex, HybridQuantity::Scalar::Ex), view(Ey, HybridQuantity::Scalar::Ey),
                       view(Ez, HybridQuantity::Scalar::Ez), view(Bx, HybridQuantity::Scalar::Bx),
                       view(By, HybridQuantity::Scalar::By), view(Bz, HybridQuantity::Scalar::Bz)});
    }

    container_t<Float> Ex, Ey, Ez, Bx, By, Bz, density, fluxX, fluxY, fluxZ;
};


template<typename Solver>
class Offloader : gpu::BaseOffloader<Solver>
{
    using This  = Offloader<Solver>;
    using Super = gpu::BaseOffloader<Solver>;

    static constexpr auto dimension       = Super::dimension;
    static constexpr auto interp_order    = Super::interp_order;
    static constexpr std::size_t TP_BLOCK = 1024; // get from dict

    using GridLayout        = typename Super::GridLayout;
    using HybridState       = typename Solver::HybridModel_t::State_t;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order, /*offload*/ true>;
    using Ions              = typename Solver::Ions;
    using ParticleArray     = typename Ions::particle_array_type;
    using Particle_t        = typename ParticleArray::value_type;
    using PartIterator      = typename ParticleArray::iterator;
    using ETC_t             = PHARE::solver::gpu_mkn::ETC<GridLayout>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher_t
        = PHARE::core::GranovPusher<dimension, PartIterator, typename ETC_t::gpu_t::EM_Interop,
                                    Interpolator, BoundaryCondition, GridLayout>;

    struct PushInfo
    {
        double timestep, mass, coef = 1;
    };

    bool static inDomain(GridLayout const& layout, Particle_t const& particle) _PHARE_FN_SIG_
    {
        return core::isIn(core::cellAsPoint(particle), layout.AMRBox());
    }

    static constexpr auto push_fn
        = [] __device__(auto& particle, PushInfo& info, auto& layout, auto& etc) {
              auto&& [density, flux, em] = etc(layout);
              assert(info.coef == 1);
              assert(info.mass > 0);
              Interpolator interpolator;
              assert(layout.AMRBox().lower[0] >= 0);
              assert(particle.iCell[0] > -1);
              Pusher_t{layout, info.timestep, info.mass}.move_in_place(
                  particle, em, interpolator, [](auto const& /*TODO*/) { return true; }, layout);
              assert(particle.iCell[0] > -6);
              if (inDomain(layout, particle))
                  interpolator.particleToMesh(particle, density, flux, layout, info.coef);
          };
    static constexpr auto push_fn_0
        = [] __device__(auto i, auto particles, auto copy, auto info, auto layout, auto etc) {
              copy[i] = particles[i]; // save state
              push_fn(particles[i], *info, *layout, *etc);
          };
    static constexpr auto push_fn_1
        = [] __device__(auto i, auto particles, auto copy, auto info, auto layout, auto etc) {
              particles[i] = copy[i]; // restore state
              push_fn(particles[i], *info, *layout, *etc);
          };
    static constexpr auto push_fn_test
        = [] __device__(auto i, auto particles, auto copy, auto info, auto layout, auto etc) {};

public:
    Offloader(PHARE::initializer::PHAREDict const& dict) {}

    auto& clear()
    {
        KLOG(TRC);
        std::apply([](auto&... v) { (v.clear(), ...); },
                   std::forward_as_tuple(etcs_, particles_, particle_copies_, patch_states,
                                         post_partition, batches));
        assert(batches.size() == 0);
        return *this;
    };

    auto& alloc(GridLayout const& layout, HybridState& hybrid_state, double timestep)
    {
        auto& patch_state = patch_states.emplace_back(layout, hybrid_state);
        for (std::size_t i = 0; i < patch_state.ions.size(); ++i)
            batches.emplace_back(_alloc_(patch_state, i, timestep))->send();
        return *this;
    };


    auto _alloc_(PatchState<GridLayout>& state, std::size_t pop_idx, double timestep)
    {
        using ParticleDevMem = kul::gpu::DeviceMem<Particle_t>;
        auto& pop            = *state.ions[pop_idx];
        particles_.emplace_back(
            std::make_shared<kul::gpu::HostMem<Particle_t>>(pop.data(), pop.size()));
        particle_copies_.emplace_back(std::make_shared<ParticleDevMem>(pop.size()));
        etcs_.emplace_back(std::make_shared<ETC_t>(state, pop_idx));
        return kul::gpu::asio::BatchMaker::make_unique(*particles_.back(), *particle_copies_.back(),
                                                       PushInfo{timestep, state.masses[pop_idx]},
                                                       state.layout, *etcs_.back());
    }

    template<bool noop = false>
    void _move0()
    {
        KLOG(TRC);
        assert(patch_states.size());

        _send_em();

        for (auto& batch : batches)
        {
            batch->sync();
            if constexpr (noop)
                kul::gpu::asio::Launcher{TP_BLOCK}.send(false).launch(push_fn_test, *batch);
            else
                kul::gpu::asio::Launcher{TP_BLOCK}.send(false).launch(push_fn_0, *batch);
        }
    }

    void _send_em()
    {
        assert(patch_states.size() == batches.size());
        for (std::size_t i = 0; i < batches.size(); ++i)
        {
            auto& etc         = std::get<3>(*batches[i]);
            auto& patch_state = patch_states[i];
            etc.Ex.send(patch_state.electromag[0]);
            etc.Ey.send(patch_state.electromag[1]);
            etc.Ez.send(patch_state.electromag[2]);
            etc.Bx.send(patch_state.electromag[3]);
            etc.By.send(patch_state.electromag[4]);
            etc.Bz.send(patch_state.electromag[5]);
        }
    }

    void _move1()
    {
        KLOG(TRC);
        assert(batches.size() > 0);

        _send_em();

        visit([&](auto& batch, auto& patch_state, auto& pop, auto pop_idx) {
            auto& etc = std::get<3>(batch);
            etc.fluxX.send(patch_state.flux[0 + (3 * pop_idx)]);
            etc.fluxY.send(patch_state.flux[1 + (3 * pop_idx)]);
            etc.fluxZ.send(patch_state.flux[2 + (3 * pop_idx)]);
            etc.density.send(patch_state.density[pop_idx]);
        });

        for (auto& batch : batches)
            kul::gpu::asio::Launcher{TP_BLOCK}.send(false).launch(push_fn_1, *batch);
    }

    template<std::size_t index>
    auto& move()
    {
        if constexpr (index == 0)
            _move0();
        else if constexpr (index == 1)
            _move1();
        else
        {
            assert(false);
        }
        return *this;
    }


    void _partition_particles()
    {
        for (auto& patch_state : patch_states)
        {
            auto& layout = patch_state.layout;
            for (auto* pop : patch_state.ions)
            {
                auto range = makeRange(*pop);
                post_partition[pop->data()]
                    = std::partition(std::begin(range), std::end(range),
                                     [&](auto const& part) { return inDomain(layout, part); });
            }
        }
    }

    template<typename Fn>
    void visit(Fn&& fn)
    {
        assert(batches.size());
        std::size_t state_idx = 0, pop_idx = 0;
        std::size_t n_pops = patch_states[0].ions.size();
        for (auto& batch : batches)
        {
            fn(*batch, patch_states[state_idx], *patch_states[state_idx].ions[pop_idx], pop_idx);
            ++pop_idx;
            if (pop_idx == n_pops)
            {
                pop_idx = 0;
                ++state_idx;
            }
        }
    }

    void _copy_back_particles()
    {
        KLOG(TRC);
        visit([&](auto& batch, auto& patch_state, auto& pop, auto pop_idx) {
            auto& copy_back = batch.get(0);
            assert(pop.size() > 0);
            assert(copy_back.size() == pop.size());
            std::copy(copy_back.data(), copy_back.data() + copy_back.size(), pop.data());
        });
    }


    void _copy_back_fields()
    {
        visit([&](auto& batch, auto& patch_state, auto& pop, auto pop_idx) {
            auto& etc = std::get<3>(batch);
            etc.fluxX.take(patch_state.flux[0 + (3 * pop_idx)]);
            etc.fluxY.take(patch_state.flux[1 + (3 * pop_idx)]);
            etc.fluxZ.take(patch_state.flux[2 + (3 * pop_idx)]);
            etc.density.take(patch_state.density[pop_idx].data());
        });
    }


    auto& join()
    {
        for (auto& batch : batches)
            batch->sync();
        _copy_back_fields();
        _copy_back_particles();
        _partition_particles();
        return *this;
    }

    auto const& new_end(ParticleArray& array)
    {
        assert(array.size());
        KLOG(INF) << array.data() << " " << array.size();
        assert(post_partition.count(array.data()));
        return post_partition.at(array.data());
    }


private:
    std::unordered_map<Particle_t const*, PartIterator> post_partition;
    std::vector<PatchState<GridLayout>> patch_states;
    std::vector<std::shared_ptr<kul::gpu::HostMem<Particle_t>>> particles_;
    std::vector<std::shared_ptr<kul::gpu::DeviceMem<Particle_t>>> particle_copies_;
    std::vector<std::shared_ptr<ETC_t>> etcs_;

    std::vector<kul::gpu::HostMem<Particle_t>> particles;

    using Batch_t = std::decay_t<std::result_of_t<decltype (&This::_alloc_)(
        This, PatchState<GridLayout>&, std::size_t, double)>>;
    std::vector<Batch_t> batches;
};
} // namespace PHARE::solver::gpu_mkn

namespace PHARE::solver::gpu
{
template<typename Solver>
using Offloader = PHARE::solver::gpu_mkn::Offloader<Solver>;
} // namespace PHARE::solver::gpu

#endif // defined(PHARE_WITH_GPU_MKN)




#if defined(PHARE_WITH_GPU_THRUST)
namespace PHARE::solver::gpu_thrust
{
template<typename Solver>
class Offloader : gpu::BaseOffloader<Solver, Offloader<Solver>>
{
    using Super                        = gpu::BaseOffloader<Solver, Offloader<Solver>>;
    static constexpr auto dimension    = Super::dimension;
    static constexpr auto interp_order = Super::interp_order;

    Offloader(Solver& solver, PHARE::initializer::PHAREDict const& dict)
        : Super{solver}
    {
    }

    template<typename T>
    void post(std::vector<T> const& v) override
    {
    }

    template<typename T>
    void post(T* const* data, std::size_t size) override
    {
    }


    template<typename T>
    T* get()
    {
    }

    template<typename T>
    T* alloc(std::vector<T> const&) override
    {
    }
    template<typename T>
    T* alloc(std::size_t) override
    {
    }

private:
};
} // namespace PHARE::solver::gpu_thrust


namespace PHARE::solver::gpu
{
template<typename Solver>
using Offloader = PHARE::solver::gpu_thrust::Offloader<Solver>;
} // namespace PHARE::solver::gpu
#endif // defined(PHARE_WITH_GPU_THRUST)



#endif /*PHARE_SOLVER_GPU_OFFLOADER_H*/
