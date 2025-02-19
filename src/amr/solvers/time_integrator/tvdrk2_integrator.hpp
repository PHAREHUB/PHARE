#ifndef PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class TVDRK2Integrator
{
    using VecFieldT = typename MHDModel::vecfield_type;
    using MHDStateT = typename MHDModel::state_type;

    using Layout        = typename MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<Layout>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

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
        tvdrk2_step_(layouts, state, RKPair_t{w0_, state}, RKPair_t{w1_, state1_});
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state1_); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state1_);
    }

private:
    static constexpr auto w0_{0.5};
    static constexpr auto w1_{0.5};

    Euler<FVMethodStrategy, MHDModel> euler_;
    RKUtils_t tvdrk2_step_;

    MHDStateT state1_{"state1"};
};

} // namespace PHARE::solver

#endif
