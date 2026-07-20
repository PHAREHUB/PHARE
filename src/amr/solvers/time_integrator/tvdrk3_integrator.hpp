#ifndef PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"
#include "amr/solvers/solver_mhd_field_evolvers.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<typename FVMethodStrategy, typename MHDModel,
         typename MessengerT = amr::MHDMessenger<MHDModel>>
class TVDRK3Integrator : public BaseMHDTimestepper<MHDModel, MessengerT>
{
    using Super = BaseMHDTimestepper<MHDModel, MessengerT>;

    using GridLayoutT = MHDModel::gridlayout_type;
    using VecFieldT   = MHDModel::vecfield_type;
    using MHDStateT   = MHDModel::state_type;

    using Dispatchers_t = Dispatchers<MHDModel>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

public:
    TVDRK3Integrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict, /*n_extra_states=*/2}
        , euler_{dict}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, Super::MHDStateT& state,
                    Super::FluxT& fluxes, Super::Messenger& bc,
                    Super::level_t& level, double const currentTime,
                    double const newTime) override
    {
        auto& state1 = this->extra_states_[0];
        auto& state2 = this->extra_states_[1];

        this->resetButcherFluxes_(model, level);

        // U1 = Euler(Un)
        euler_(model, state, state1, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state.E, fluxes, level, w01_ * w11_);

        // U1 = Euler(U1)
        euler_(model, state1, state1, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state1.E, fluxes, level, w01_ * w11_);

        // U2 = 0.75*Un + 0.25*U1
        RKUtils_t{level, model}(newTime, state2, RKPair_t{w00_, state}, RKPair_t{w01_, state1});

        // U2 = Euler(U2)
        euler_(model, state2, state2, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state2.E, fluxes, level, w11_);

        euler_using_butcher_fluxes_(model, state, state, this->butcherE_, this->butcherFluxes_, bc,
                                    level, newTime, newTime - currentTime);
    }

    void registerResources(MHDModel& model) override
    {
        Super::registerResources(model);
        euler_.registerResources(model);
    }

    void allocate(MHDModel& model, SAMRAI::hier::Patch& patch,
                  double const allocateTime) const override
    {
        Super::allocate(model, patch, allocateTime);
        euler_.allocate(model, patch, allocateTime);
    }

    using Super::exposeFluxes;

private:
    static constexpr auto w00_{0.75};
    static constexpr auto w01_{0.25};
    static constexpr auto w10_{1. / 3.};
    static constexpr auto w11_{2. / 3.};

    Euler<FVMethodStrategy, MHDModel> euler_;
    EulerUsingComputedFlux<MHDModel> euler_using_butcher_fluxes_;
};

} // namespace PHARE::solver

#endif
