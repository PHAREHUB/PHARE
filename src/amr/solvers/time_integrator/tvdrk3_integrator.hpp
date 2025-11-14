#ifndef PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class TVDRK3Integrator : public BaseMHDTimestepper<MHDModel>
{
    using Super = BaseMHDTimestepper<MHDModel>;

    using level_t     = typename MHDModel::level_t;
    using FieldT      = typename MHDModel::field_type;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayoutT = typename MHDModel::gridlayout_type;
    using MHDStateT   = typename MHDModel::state_type;

    using Dispatchers_t = Dispatchers<GridLayoutT>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

public:
    TVDRK3Integrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
        , euler_{dict}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, level_t& level,
                    double const currentTime, double const newTime)
    {
        this->resetButcherFluxes_(model, level);

        // U1 = Euler(Un)
        euler_(model, state, state1_, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state.E, fluxes, level, w01_ * w11_);

        // U1 = Euler(U1)
        euler_(model, state1_, state1_, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state1_.E, fluxes, level, w01_ * w11_);

        // U2 = 0.75*Un + 0.25*U1
        tvdrk3_step_(level, model, newTime, state2_, RKPair_t{w00_, state},
                     RKPair_t{w01_, state1_});

        // U2 = Euler(U2)
        euler_(model, state2_, state2_, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state2_.E, fluxes, level, w11_);

        euler_using_butcher_fluxes_(model, state, state, this->butcherE_, this->butcherFluxes_, bc,
                                    level, newTime, newTime - currentTime);

        // Un+1 = 1/3*Un + 2/3*Euler(U2)
        // tvdrk3_step_(level, model, newTime, state, RKPair_t{w10_, state}, RKPair_t{w11_,
        // state2_});
    }

    void registerResources(MHDModel& model)
    {
        Super::registerResources(model);
        model.resourcesManager->registerResources(state1_);
        model.resourcesManager->registerResources(state2_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        Super::allocate(model, patch, allocateTime);
        model.resourcesManager->allocate(state1_, patch, allocateTime);
        model.resourcesManager->allocate(state2_, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const
    {
        auto fill_info = [&](auto& state) {
            info.ghostDensity.push_back(state.rho.name());
            info.ghostMomentum.push_back(state.rhoV.name());
            info.ghostTotalEnergy.push_back(state.Etot.name());
            info.ghostElectric.push_back(state.E.name());
            info.ghostMagnetic.push_back(state.B.name());
            info.ghostCurrent.push_back(state.J.name());
        };

        fill_info(state1_);
        fill_info(state2_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::tuple_cat(Super::getCompileTimeResourcesViewList(),
                              std::forward_as_tuple(state1_, state2_));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::tuple_cat(Super::getCompileTimeResourcesViewList(),
                              std::forward_as_tuple(state1_, state2_));
    }

    using Super::exposeFluxes;

private:
    static constexpr auto w00_{0.75};
    static constexpr auto w01_{0.25};
    static constexpr auto w10_{1. / 3.};
    static constexpr auto w11_{2. / 3.};

    Euler<FVMethodStrategy, MHDModel> euler_;
    EulerUsingComputedFlux<MHDModel> euler_using_butcher_fluxes_;
    RKUtils_t tvdrk3_step_;

    MHDStateT state1_{"state1"};
    MHDStateT state2_{"state2"};
};

} // namespace PHARE::solver

#endif
