#ifndef PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "core/numerics/time_integrator/euler.hpp"
#include "core/models/mhd_state.hpp"
#include "core/numerics/time_integrator/time_integrator_utils.hpp"

namespace PHARE::core
{
template<template<typename> typename FVMethod, typename MHDModel>
class TVDRK3Integrator
{
    using VecFieldT = typename MHDModel::vecfield_type;
    using MHDStateT = typename MHDModel::state_type;

    using ModelView_t = solver::MHDModelView<MHDModel>;
    using RKUtils_t   = ModelView_t::RKUtils_t;

public:
    TVDRK3Integrator(PHARE::initializer::PHAREDict const& dict)
        : euler_{dict}
    {
    }

    void operator()(auto& layouts, auto& state, auto& fluxes, auto& bc, auto& level,
                    double const currentTime, double const newTime)
    {
        double const dt = newTime - currentTime;

        // U1 = Euler(Un)
        euler_(layouts, state, state1_, fluxes, bc, level, currentTime, newTime);

        // U1 = Euler(U1)
        euler_(layouts, state1_, state1_, fluxes, bc, level, currentTime, newTime);

        // U2 = 0.75*Un + 0.25*U1
        tvdrk3_step_(layouts, state2_, std::make_pair(w00_, state), std::make_pair(w01_, state1_));

        // U2 = Euler(U2)
        euler_(layouts, state2_, state2_, fluxes, bc, level, currentTime, newTime);

        // Un+1 = 1/3*Un + 2/3*Euler(U2)
        tvdrk3_step_(layouts, state, std::make_pair(w10_, state), std::make_pair(w11_, state2_));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(state1_, state2_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state1_, state2_);
    }

private:
    double const w00_{0.75};
    double const w01_{0.25};
    double const w10_{1. / 3.};
    double const w11_{2. / 3.};

    Euler<FVMethod, MHDModel> euler_;
    RKUtils_t tvdrk3_step_;

    MHDStateT state1_{"state1"};
    MHDStateT state2_{"state2"};
};

} // namespace PHARE::core

#endif
