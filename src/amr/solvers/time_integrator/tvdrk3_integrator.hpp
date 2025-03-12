#ifndef PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP

#include "core/data/vecfield/vecfield.hpp"
#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class TVDRK3Integrator
{
    using VecFieldT = typename MHDModel::vecfield_type;
    using MHDStateT = typename MHDModel::state_type;

    using Layout        = typename MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<Layout>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

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
        tvdrk3_step_(layouts, state2_, RKPair_t{w00_, state}, RKPair_t{w01_, state1_});

        // U2 = Euler(U2)
        euler_(layouts, state2_, state2_, fluxes, bc, level, currentTime, newTime);

        // Un+1 = 1/3*Un + 2/3*Euler(U2)
        tvdrk3_step_(layouts, state, RKPair_t{w10_, state}, RKPair_t{w11_, state2_});
    }

    void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(state1_);
        model.resourcesManager->registerResources(state2_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(state1_, patch, allocateTime);
        model.resourcesManager->allocate(state2_, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const
    {
        info.ghostDensity.push_back(state1_.rho.name());
        info.ghostVelocity.push_back(core::VecFieldNames{state1_.V});
        info.ghostPressure.push_back(state1_.P.name());
        info.ghostElectric.push_back(core::VecFieldNames{state1_.E});

        info.ghostDensity.push_back(state2_.rho.name());
        info.ghostVelocity.push_back(core::VecFieldNames{state2_.V});
        info.ghostPressure.push_back(state2_.P.name());
        info.ghostElectric.push_back(core::VecFieldNames{state2_.E});
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
    static constexpr auto w00_{0.75};
    static constexpr auto w01_{0.25};
    static constexpr auto w10_{1. / 3.};
    static constexpr auto w11_{2. / 3.};

    Euler<FVMethodStrategy, MHDModel> euler_;
    RKUtils_t tvdrk3_step_;

    MHDStateT state1_{"state1"};
    MHDStateT state2_{"state2"};
};

} // namespace PHARE::solver

#endif
