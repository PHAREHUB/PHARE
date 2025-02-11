#ifndef PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "core/models/mhd_state.hpp"
#include "core/numerics/time_integrator/euler.hpp"
#include "core/numerics/time_integrator/time_integrator_utils.hpp"
#include <utility>

namespace PHARE::core
{
template<template<typename> typename FVMethod, typename MHDModel>
class TVDRK2Integrator
{
    using VecFieldT = typename MHDModel::vecfield_type;
    using MHDStateT = typename MHDModel::state_type;

    using ModelView_t = solver::MHDModelView<MHDModel>;
    using RKUtils_t   = ModelView_t::RKUtils_t;

public:
    TVDRK2Integrator(PHARE::initializer::PHAREDict const& dict)
        : euler_{dict}
    {
    }

    void operator()(auto& layouts, MHDStateT& state, auto& fluxes, auto& bc, auto& level,
                    double const currentTime, double const newTime)
    {
        double const dt = newTime - currentTime;

        // U1 = Euler(Un)
        euler_(layouts, state, state1_, fluxes, bc, level, currentTime, newTime);

        // U1 = Euler(U1)
        euler_(layouts, state1_, state1_, fluxes, bc, level, currentTime, newTime);

        // Un+1 = 0.5*Un + 0.5*Euler(U1)
        tvdrk2_step_(layouts, state, std::make_pair(w0_, state), std::make_pair(w1_, state1_));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state1_); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state1_);
    }

private:
    double const w0_{0.5};
    double const w1_{0.5};

    Euler<FVMethod, MHDModel> euler_;
    RKUtils_t tvdrk2_step_;

    MHDStateT state1_{"state1"};
};

} // namespace PHARE::core

#endif
