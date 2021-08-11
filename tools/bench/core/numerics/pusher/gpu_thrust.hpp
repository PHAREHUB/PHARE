#ifndef PHARE_GPU_THRUST_HPP
#define PHARE_GPU_THRUST_HPP

#include "kul/log.hpp"

#include "benchmark/benchmark.h"


#include <thrust/generate.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/host_vector.h>

#include <algorithm>
#include <cstdlib>


#include "core/numerics/pusher/boris.h"
#include "core/numerics/pusher/granov.h"


#include "tests/simulator/per_test.h"

namespace PHARE::gpu_thrust
{
template<typename PHARE_TYPES>
struct push_functor : public thrust::unary_function<PHARE::core::Particle<PHARE_TYPES::dimension>,
                                                    PHARE::core::Particle<PHARE_TYPES::dimension>>
{
    static constexpr auto dim    = PHARE_TYPES::dimension;
    static constexpr auto interp = PHARE_TYPES::interp_order;

    using GridLayout_t = typename PHARE_TYPES::GridLayout_t;

    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Ions_t            = typename PHARE_TYPES::Ions_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using Particle          = typename ParticleArray::value_type;
    using PartIterator      = typename ParticleArray::iterator;

    using Float = double;


    struct __F__
    {
        __F__(thrust::device_vector<Float>&& f0)
            : f{f0}
            , p{f0.data()}
        {
        }
        auto __device__ operator()(std::uint32_t i) { return p[i]; }
        auto __device__ operator()(std::uint32_t i) const { return p[i]; }

        thrust::device_vector<Float> f;
        thrust::device_ptr<Float> p;
    };

    struct __VF__
    {
        auto __device__ getComponents() { return std::forward_as_tuple(x, y, z); }
        auto __device__ getComponents() const { return std::forward_as_tuple(x, y, z); }

        __F__ x, y, z;
    };
    struct __EM__
    {
        __VF__ E, B;
    };

    using Electromag_t = __EM__;
    using Pusher_t     = PHARE::core::GranovPusher<dim, PartIterator, Electromag_t, Interpolator,
                                               BoundaryCondition, GridLayout_t>;

    template<typename EM_input>
    push_functor(GridLayout_t layout_, EM_input const& EM, double timestep_)
        : timestep{timestep_}
        , layout(layout_)
        , em{{{std::vector<double>(EM[0].data(), EM[0].data() + EM[0].size())},
              {std::vector<double>(EM[1].data(), EM[1].data() + EM[1].size())},
              {std::vector<double>(EM[2].data(), EM[2].data() + EM[2].size())}},
             {{std::vector<double>(EM[3].data(), EM[3].data() + EM[3].size())},
              {std::vector<double>(EM[4].data(), EM[4].data() + EM[4].size())},
              {std::vector<double>(EM[5].data(), EM[5].data() + EM[5].size())}}}
    {
    }

    auto operator()(PHARE::core::Particle<dim>& particle) _PHARE_FN_SIG_
    {
        Interpolator interpolator;
        Pusher_t pusher{layout, timestep, /*mass=*/1};
        pusher.move_in_place(
            /*Particle_t&*/ particle,                  //
            /*Electromag const&*/ em,                  //
            /*Interpolator&*/ interpolator,            //
            [](auto const& /*part*/) { return true; }, //
            /*GridLayout const&*/ layout);
        return particle;
    }


    double timestep = 0;
    GridLayout_t layout;
    __EM__ em;

    template<typename GPU_setup>
    static void push(benchmark::State& benchmark_state, double timestep)
    {
        std::string job_id = "job_" + std::to_string(1) + "d";

        for (auto _ : benchmark_state)
        {
            GPU_setup setup{job_id};
            KLOG(INF) << "GPU PARTICLES: " << setup.packer.n_particles;
            auto& state = setup.states[0];
            thrust::transform(state.ions[0]->begin(), state.ions[0]->end(), state.ions[0]->begin(),
                              push_functor{state.layout, state.electromag, timestep});
            KLOG(INF) << " GPU particle  0 delta after  " << state.ions[0]->front().delta[0];
        }
    }
};

} // namespace PHARE::gpu_thrust

#endif /*PHARE_GPU_THRUST_HPP*/
