

#include "range.hpp"
#include "kul/log.hpp"
#include "pusher_bench.h"

#include "gpu_thrust.hpp"

#include <algorithm>

#include "benchmark/benchmark.h"

#include "core/numerics/pusher/boris.h"
#include "core/numerics/pusher/granov.h"

#include "tools/hw/gpu.hpp"

using namespace PHARE::core::bench;

static constexpr std::size_t X = 1024, Y = 1024, Z = 10;        // 40;
static constexpr std::size_t TPB_X = 16, TPB_Y = 16, TPB_Z = 1; // 4;
static constexpr std::size_t MAX_PARTICLES = X * Y * Z;
static constexpr std::uint32_t TP_BLOCK    = 1024;

static constexpr std::size_t PUSH_TIMES = 2;

static constexpr double TIMESTEP = .0000001; // TORM (somehow)


template<typename PHARE_TYPES>
class Setup
{
    using Hierarchy_t   = typename PHARE_TYPES::hierarchy_t;
    using HybridModel_t = typename PHARE_TYPES::HybridModel_t;
    using HybridState_t = typename HybridModel_t::State_t;
    using Types         = PHARE_TYPES;

    static constexpr auto dim           = PHARE_TYPES::dimension;
    static constexpr auto interp        = PHARE_TYPES::interp_order;
    static auto constexpr nbRefineParts = PHARE_TYPES::nbRefinedPart;

    auto static fill_states(Setup& self)
    {
        std::vector<PHARE::gpu::PatchState<PHARE_TYPES>> states;
        PHARE::amr::visitHierarchy<typename PHARE_TYPES::GridLayout_t>(
            self.hierarchy, *self.hybridModel.resourcesManager,
            [&](auto& gridLayout, std::string patchID, size_t) {
                states.emplace_back(gridLayout, self.state, patchID);
            },
            self.topLvl, self.topLvl + 1, self.hybridModel);
        return states;
    }

public:
    Setup(std::string job_id)
        : sim{job_id}
    {
    }

    SimulatorTestParam<dim, interp, nbRefineParts> sim;
    Hierarchy_t& hierarchy{*sim.hierarchy};
    HybridModel_t& hybridModel{*sim.getHybridModel()};
    HybridState_t& state{hybridModel.state};
    int topLvl{hierarchy.getNumberOfLevels() - 1};

    std::vector<PHARE::gpu::PatchState<PHARE_TYPES>> states{fill_states(*this)};
};

template<typename PHARE_TYPES, bool async = false>
class GPU_setup : public Setup<PHARE_TYPES>
{
public:
    using Super = Setup<PHARE_TYPES>;
    using Super::hierarchy;
    using Super::hybridModel;
    using Super::state;
    using Super::states;
    using Super::topLvl;

    GPU_setup(std::string job_id)
        : Super{job_id}
    {
    }

    PHARE::gpu::ParticlePatchState<PHARE_TYPES> packer{states, async};
};


template<typename PHARE_TYPES>
auto per_particle_alloc()
{
    static constexpr auto dim = PHARE_TYPES::dimension;

    return (MAX_PARTICLES * sizeof(PHARE::gpu::PatchStatePerParticle<PHARE_TYPES, true>)
            + MAX_PARTICLES * sizeof(PHARE::core::Particle<dim>));
}


template<typename PHARE_TYPES>
struct Pusher
{
    static constexpr auto dim    = PHARE_TYPES::dimension;
    static constexpr auto interp = PHARE_TYPES::interp_order;

    using Interpolator      = PHARE::core::Interpolator<dim, interp, /*offload*/ true>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Ions_t            = typename PHARE_TYPES::Ions_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using Particle          = typename ParticleArray::value_type;
    using PartIterator      = typename ParticleArray::iterator;

    struct CPU
    {
        using Electromag_t = typename PHARE_TYPES::Electromag_t;
        using GridLayout_t = typename PHARE_TYPES::GridLayout_t;
        using BorisPusher_t
            = PHARE::core::BorisPusher<dim, PartIterator, Electromag_t, Interpolator,
                                       BoundaryCondition, GridLayout_t>;

        static void push(benchmark::State& benchmark_state)
        {
            std::string job_id = "job_" + std::to_string(1) + "d";

            KLOG(INF) << "MAX_PARTICLES spawned: " << MAX_PARTICLES;
            KLOG(INF) << "ALL_PARTICLES bytes  : " << per_particle_alloc<PHARE_TYPES>();

            for (auto _ : benchmark_state)
            {
                Setup<PHARE_TYPES> setup{job_id};
                Interpolator interpolator;


                PHARE::amr::visitHierarchy<GridLayout_t>(
                    setup.hierarchy, *setup.hybridModel.resourcesManager,
                    [&](auto& layout, std::string patchID, size_t) {
                        for (auto& pop : setup.state.ions)
                        {
                            auto range = PHARE::core::makeRange(pop.domainParticles());
                            BorisPusher_t pusher;
                            pusher.setMeshAndTimeStep(layout.meshSize(), TIMESTEP);
                            for (std::size_t j = 0; j < PUSH_TIMES; ++j)
                                pusher.move(
                                    range, range, setup.state.electromag, 1, interpolator,
                                    [](auto const& /*part*/) { return true; }, layout);
                        }
                    },
                    setup.topLvl, setup.topLvl + 1, setup.hybridModel);


                // if (PHARE::core::mpi::any(PHARE::core::Errors::instance().any()))
                //     throw std::runtime_error("errors");
            }
        }

        template<typename Setup>
        static void push(Setup& setup, PHARE::gpu::PatchState<PHARE_TYPES>& patch_state,
                         bool first = true)
        {
            Interpolator interpolator;

            BorisPusher_t pusher;
            pusher.setMeshAndTimeStep(patch_state.layout.meshSize(), TIMESTEP);

            PHARE::amr::visitHierarchy<GridLayout_t>(
                setup.hierarchy, *setup.hybridModel.resourcesManager,
                [&](auto& layout, std::string patchID, size_t) {
                    if (patchID == patch_state.id)
                    {
                        auto& pop = *[&](auto& ions) {
                            if (first)
                                return ions.front();
                            return ions.back();
                        }(patch_state.ions);

                        if (first)
                            KLOG(INF) << "CPU particle  0 delta before  " << pop.front().delta[0];
                        else
                            KLOG(INF) << "CPU particle -1 delta before  " << pop.back().delta[0];

                        auto range = PHARE::core::makeRange(pop);
                        for (std::size_t j = 0; j < PUSH_TIMES; ++j)
                            pusher.move(
                                range, range, setup.state.electromag, 1, interpolator,
                                [](auto const& /*part*/) { return true; }, layout);

                        if (first)
                            KLOG(INF) << "CPU particle  0 delta after  " << pop.front().delta[0];
                        else
                            KLOG(INF) << "CPU particle -1 delta after  " << pop.back().delta[0];
                    }
                },
                setup.topLvl, setup.topLvl + 1, setup.hybridModel);

            if (PHARE::core::mpi::any(PHARE::core::Errors::instance().any()))
                throw std::runtime_error("errors");
        }

        template<typename Setup>
        static void push(Setup& setup)
        {
            push(setup, setup.states.front());
            if (setup.states.size() > 1)
                push(setup, setup.states.back(), false);
        }
    };

    struct Sync
    {
        static void gpu_push(PHARE::gpu::PatchStatePerParticle<PHARE_TYPES, true>* ppsp) __global__
        {
            auto i = kul::gpu::idx();

            if (i >= ppsp->n_particles())
                return;

            auto patchStateIDX = (*ppsp)[i];
            auto& layout       = ppsp->gridLayouts->layouts[patchStateIDX];
            auto electromag    = ppsp->electromags->electromag(patchStateIDX);

            using Electromag_t = decltype(electromag);
            using GridLayout_t = std::decay_t<decltype(layout)>;
            using Pusher_t
                = PHARE::core::GranovPusher<dim, PartIterator, Electromag_t, Interpolator,
                                            BoundaryCondition, GridLayout_t>;

            Interpolator interpolator;
            Pusher_t pusher{layout, TIMESTEP, /*mass=*/1};
            for (std::size_t j = 0; j < PUSH_TIMES; ++j)
                pusher.move_in_place(
                    /*Particle_t&*/ ppsp->particles->particles[i],                         //
                    /*Electromag const&*/ electromag,                                      //
                    /*Interpolator&*/ interpolator,                                        //
                    /*ParticleSelector const&*/ [](auto const& /*part*/) { return true; }, //
                    /*GridLayout const&*/ layout                                           //
                );
        }

        static void push(benchmark::State& benchmark_state)
        {
            std::string job_id = "job_" + std::to_string(1) + "d";

            KLOG(INF) << "MAX_PARTICLES spawned: " << MAX_PARTICLES;
            KLOG(INF) << "ALL_PARTICLES bytes  : " << per_particle_alloc<PHARE_TYPES>();

            for (auto _ : benchmark_state)
            {
                GPU_setup<PHARE_TYPES> setup{job_id};
                KLOG(INF) << "GPU PARTICLES: " << setup.packer.n_particles;

                kul::gpu::Launcher{X, Y, Z, TPB_X, TPB_Y, TPB_Z}(gpu_push, setup.packer());

                auto gpu_particles = setup.packer.particles->particles();
                KLOG(INF) << " GPU particle  0 delta after  " << gpu_particles.front().delta[0];
                KLOG(INF) << " GPU particle -1 delta after  " << gpu_particles.back().delta[0];

                CPU::push(setup);
            }
        }
    };


    struct Async
    {
        using Electromag_t = PHARE::gpu::Electromag<PHARE_TYPES, /*gpu=*/true>;
        using GridLayout_t = typename PHARE_TYPES::GridLayout_t;
        using Pusher_t
            = PHARE::core::GranovPusher<dim, PartIterator, typename Electromag_t::Interop,
                                        Interpolator, BoundaryCondition, GridLayout_t>;

        static void push(benchmark::State& benchmark_state)
        {
            std::string job_id = "job_" + std::to_string(1) + "d";

            for (auto _ : benchmark_state)
            {
                GPU_setup<PHARE_TYPES, /*async=*/true> setup{job_id};
                KLOG(INF) << "GPU PARTICLES: " << setup.packer.n_particles;

                auto& states          = setup.states;
                std::size_t state_idx = 0, pop_idx = 0;
                std::size_t n_pops = states[0].ions.size();

                auto batches = PHARE::generate(
                    [&](auto&&) {
                        auto& state = states[state_idx];
                        auto& pop   = *state.ions[pop_idx];
                        ++pop_idx;

                        if (pop_idx == n_pops)
                        {
                            pop_idx = 0;
                            ++state_idx;
                        }
                        return kul::gpu::asio::Launcher{TP_BLOCK, 1}(
                            [] __device__(auto i, auto particles, auto layout, auto em) {
                                Interpolator interpolator;
                                Pusher_t pusher{*layout, TIMESTEP, /*mass=*/1};
//                                 for (std::size_t j = 0; j < PUSH_TIMES; ++j)
//                                     pusher.move_in_place(
//                                         /*Particle_t&*/ particles[i],              //
//                                         /*Electromag const&*/ (*em)(),             //
//                                         /*Interpolator&*/ interpolator,            //
//                                         [](auto const& /*part*/) { return true; }, //
//                                         /*GridLayout const&*/ *layout);
                            },
                            pop, state.layout, PHARE::gpu::Electromag<PHARE_TYPES>{state});
                    },
                    states.size() * n_pops)();

                KLOG(INF) << " GPU particle  0 delta after  "
                          << batches.front()->get(0).front().delta[0];
                KLOG(INF) << " GPU particle -1 delta after  "
                          << batches.back()->get(0).back().delta[0];
            }
        }
    };

    struct AsyncPinned
    {
        using Electromag_t = PHARE::gpu::Electromag<PHARE_TYPES, /*gpu=*/true>;
        using GridLayout_t = typename PHARE_TYPES::GridLayout_t;
        using Pusher_t
            = PHARE::core::GranovPusher<dim, PartIterator, typename Electromag_t::Interop,
                                        Interpolator, BoundaryCondition, GridLayout_t>;

        static void push(benchmark::State& benchmark_state, std::size_t n_batches = 1)
        {
            std::string job_id = "job_" + std::to_string(1) + "d";

            for (auto _ : benchmark_state)
            {
                GPU_setup<PHARE_TYPES, /*async=*/true> setup{job_id};
                KLOG(INF) << "GPU PARTICLES: " << setup.packer.n_particles;

                auto& states          = setup.states;
                std::size_t state_idx = 0, pop_idx = 0;
                std::size_t n_pops = states[0].ions.size();
                std::vector<kul::gpu::HostMem<Particle>> pinned;
                auto batches = PHARE::generate(
                    [&](auto&&) {
                        auto& state = states[state_idx];
                        auto& pop   = *state.ions[pop_idx];
                        ++pop_idx;

                        if (pop_idx == n_pops)
                        {
                            pop_idx = 0;
                            ++state_idx;
                        }
                        return kul::gpu::asio::Launcher{TP_BLOCK, n_batches}(
                            [] __device__(auto i, auto particles, auto layout, auto em) {
//                                 Interpolator interpolator;
//                                 Pusher_t pusher{*layout, TIMESTEP, /*mass=*/1};
//                                 for (std::size_t j = 0; j < PUSH_TIMES; ++j)
//                                     pusher.move_in_place(
//                                         /*Particle_t&*/ particles[i],              //
//                                         /*Electromag const&*/ (*em)(),             //
//                                         /*Interpolator&*/ interpolator,            //
//                                         [](auto const& /*part*/) { return true; }, //
//                                         /*GridLayout const&*/ *layout);
                            },
                            pinned.emplace_back(pop.data(), pop.size()), //
                            state.layout,                              //
                            PHARE::gpu::Electromag<PHARE_TYPES>{state} //
                        );
                    },
                    states.size() * n_pops)();

                KLOG(INF) << " GPU particle  0 delta after  "
                          << batches.front()->get(0).front().delta[0];
                KLOG(INF) << " GPU particle -1 delta after  "
                          << batches.back()->get(n_batches - 1).back().delta[0];
            }
        }
    };



    struct push_functor : public PHARE::gpu_thrust::push_functor<PHARE_TYPES>
    {
        using Super = PHARE::gpu_thrust::push_functor<PHARE_TYPES>;

        static void push(benchmark::State& benchmark_state)
        {
            Super::push(benchmark_state, TIMESTEP);
        }
    };
};


// first one is cold cache so ignore
template<std::size_t dim, std::size_t interp, std::size_t nbRefineParts = 2>
void cpu_push(benchmark::State& benchmark_state)
{
    Pusher<PHARE::PHARE_Types<dim, interp, nbRefineParts>>::CPU::push(benchmark_state);
}
BENCHMARK_TEMPLATE(cpu_push, /*dim=*/1, /*interp=*/3);


template<std::size_t dim, std::size_t interp, std::size_t nbRefineParts = 2>
void push_sync(benchmark::State& benchmark_state)
{
    Pusher<PHARE::PHARE_Types<dim, interp, nbRefineParts>>::Sync::push(benchmark_state);
}
BENCHMARK_TEMPLATE(push_sync, /*dim=*/1, /*interp=*/3);


// template<std::size_t dim, std::size_t interp, std::size_t nbRefineParts = 2>
// void push_async(benchmark::State& benchmark_state)
// {
//     Pusher<PHARE::PHARE_Types<dim, interp, nbRefineParts>>::Async::push(benchmark_state);
// }
// BENCHMARK_TEMPLATE(push_async, /*dim=*/1, /*interp=*/3);


template<std::size_t dim, std::size_t interp, std::size_t batches, std::size_t nbRefineParts = 2>
void push_async_pinned(benchmark::State& benchmark_state)
{
    Pusher<PHARE::PHARE_Types<dim, interp, nbRefineParts>>::AsyncPinned::push(benchmark_state,
                                                                              batches);
}
BENCHMARK_TEMPLATE(push_async_pinned, /*dim=*/1, /*interp=*/3, /*batches=*/1);
// BENCHMARK_TEMPLATE(push_async_pinned, /*dim=*/1, /*interp=*/3, /*batches=*/2);
// BENCHMARK_TEMPLATE(push_async_pinned, /*dim=*/1, /*interp=*/3, /*batches=*/4);
// BENCHMARK_TEMPLATE(push_async_pinned, /*dim=*/1, /*interp=*/3, /*batches=*/5);
// BENCHMARK_TEMPLATE(push_async_pinned, /*dim=*/1, /*interp=*/3, /*batches=*/10);

// hot cache so real
// BENCHMARK_TEMPLATE(cpu_push, /*dim=*/1, /*interp=*/3);


/*
void thrust_sort(benchmark::State& benchmark_state)
{
    for (auto _ : benchmark_state)
    { // generate 32M random numbers serially
        thrust::host_vector<int> h_vec(32 << 20);
        std::generate(h_vec.begin(), h_vec.end(), rand);

        // transfer data to the device
        thrust::device_vector<int> d_vec = h_vec;

        // sort data on the device (846M keys per second on GeForce GTX 480)
        thrust::sort(d_vec.begin(), d_vec.end());

        // transfer data back to host
        thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
    }
}
BENCHMARK(thrust_sort);

void thrust_eg(benchmark::State& benchmark_state)
{
    for (auto _ : benchmark_state)
    { // generate 32M random numbers serially
        thrust::host_vector<int> h_vec(32 << 20);
        std::generate(h_vec.begin(), h_vec.end(), rand);

        // transfer data to the device
        thrust::device_vector<int> d_vec = h_vec;

        // sort data on the device (846M keys per second on GeForce GTX 480)

        thrust::transform(d_vec.begin(), d_vec.end(), d_vec.begin(),
                          [] __device__(auto& v0) { return v0 + 1; });

        // thrust::sort(d_vec.begin(), d_vec.end());

        // transfer data back to host
        thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());
    }
}
BENCHMARK(thrust_eg);*/



template<std::size_t dim, std::size_t interp, std::size_t nbRefineParts = 2>
void thrust_push(benchmark::State& benchmark_state)
{
    using PHARE_Types = PHARE::PHARE_Types<dim, interp, nbRefineParts>;
    using GPU_SETUP   = GPU_setup<PHARE_Types, /*async=*/true>;
    Pusher<PHARE_Types>::push_functor::template push<GPU_SETUP>(benchmark_state);
}
// BENCHMARK_TEMPLATE(thrust_push, /*dim=*/1, /*interp=*/3); // takes ages?





int main(int argc, char** argv, char** envp)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
