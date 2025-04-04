#ifndef PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP

#include "core/data/vecfield/vecfield.hpp"
#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class TVDRK2Integrator
{
    using level_t   = typename MHDModel::level_t;
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

    void operator()(MHDModel& model, MHDStateT& state, auto& fluxes, auto& bc, level_t& level,
                    double const currentTime, double const newTime)
    {
        double const dt = newTime - currentTime;

        // U1 = Euler(Un)
        euler_(model, state, state1_, fluxes, bc, level, currentTime, newTime);

        // U1 = Euler(U1)
        euler_(model, state1_, state1_, fluxes, bc, level, currentTime, newTime);

        // Un+1 = 0.5*Un + 0.5*Euler(U1)
        tvdrk2_step_(level, model, state, RKPair_t{w0_, state}, RKPair_t{w1_, state1_});
    }

    void registerResources(MHDModel& model) { model.resourcesManager->registerResources(state1_); }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(state1_, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const
    {
        info.ghostDensity.push_back(state1_.rho.name());
        info.ghostVelocity.push_back(core::VecFieldNames{state1_.V});
        info.ghostPressure.push_back(state1_.P.name());
        info.ghostElectric.push_back(core::VecFieldNames{state1_.E});
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
