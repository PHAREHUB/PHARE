#ifndef PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK2_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"
#include "amr/solvers/time_integrator/euler.hpp"

namespace PHARE::solver
{
template<typename FVMethodStrategy, typename MHDModel,
         typename MessengerT = amr::MHDMessenger<MHDModel>>
class TVDRK2Integrator : public BaseMHDTimestepper<MHDModel, MessengerT>
{
    using Super = BaseMHDTimestepper<MHDModel, MessengerT>;

public:
    TVDRK2Integrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict, /*n_extra_states=*/1}
        , euler_{dict}
    {
    }

    void operator()(MHDModel& model, Super::MHDStateT& state, Super::FluxT& fluxes,
                    Super::Messenger& bc, Super::level_t& level, double const currentTime,
                    double const newTime) override
    {
        auto& state1 = this->extra_states_[0];

        this->resetButcherFluxes_(model, level);

        // U1 = Euler(Un)
        euler_(model, state, state1, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state.E, fluxes, level, w0_);

        // U1 = Euler(U1)
        euler_(model, state1, state1, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state1.E, fluxes, level, w1_);

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
    static constexpr auto w0_{0.5};
    static constexpr auto w1_{0.5};

    Euler<FVMethodStrategy, MHDModel> euler_;
    EulerUsingComputedFlux<MHDModel> euler_using_butcher_fluxes_;
};

} // namespace PHARE::solver

#endif
