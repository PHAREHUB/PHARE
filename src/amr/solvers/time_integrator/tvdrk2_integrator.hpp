#ifndef PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class TVDRK2Integrator : public BaseMHDTimestepper<MHDModel>
{
    using Super = BaseMHDTimestepper<MHDModel>;

    using level_t     = typename MHDModel::level_t;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayoutT = typename MHDModel::gridlayout_type;
    using MHDStateT   = typename MHDModel::state_type;

    using Dispatchers_t = Dispatchers<GridLayoutT>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

public:
    TVDRK2Integrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
        , euler_{dict}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, MHDStateT& state, auto& fluxes, auto& bc, level_t& level,
                    double const currentTime, double const newTime)
    {
        this->resetButcherFluxes_(model, level);

        // U1 = Euler(Un)
        euler_(model, state, state1_, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state.E, fluxes, level, w1_);

        // U1 = Euler(U1)
        euler_(model, state1_, state1_, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state1_.E, fluxes, level, w1_);

        euler_using_butcher_fluxes_(model, state, state, this->butcherE_, this->butcherFluxes_, bc,
                                    level, newTime, newTime - currentTime);

        // Un+1 = 0.5*Un + 0.5*Euler(U1)
        // tvdrk2_step_(level, model, newTime, state, RKPair_t{w0_, state}, RKPair_t{w1_, state1_});
    }

    void registerResources(MHDModel& model)
    {
        Super::registerResources(model);
        model.resourcesManager->registerResources(state1_);
        euler_.registerResources(model);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        Super::allocate(model, patch, allocateTime);
        model.resourcesManager->allocate(state1_, patch, allocateTime);
        euler_.allocate(model, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const
    {
        info.ghostDensity.push_back(state1_.rho.name());
        info.ghostMomentum.push_back(state1_.rhoV.name());
        info.ghostTotalEnergy.push_back(state1_.Etot.name());
        info.ghostElectric.push_back(state1_.E.name());
        info.ghostMagnetic.push_back(state1_.B.name());
        info.ghostCurrent.push_back(state1_.J.name());
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::tuple_cat(Super::getCompileTimeResourcesViewList(),
                              std::forward_as_tuple(state1_));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::tuple_cat(Super::getCompileTimeResourcesViewList(),
                              std::forward_as_tuple(state1_));
    }

    using Super::exposeFluxes;

private:
    static constexpr auto w0_{0.5};
    static constexpr auto w1_{0.5};

    Euler<FVMethodStrategy, MHDModel> euler_;
    EulerUsingComputedFlux<MHDModel> euler_using_butcher_fluxes_;
    RKUtils_t tvdrk2_step_;

    MHDStateT state1_{"state1"};
};

} // namespace PHARE::solver

#endif
