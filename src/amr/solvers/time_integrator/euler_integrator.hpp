#ifndef PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"

#include "amr/solvers/time_integrator/euler.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"

namespace PHARE::solver
{


template<typename FVMethodStrategy, typename MHDModel,
         typename MessengerT = amr::MHDMessenger<MHDModel>>
class EulerIntegrator : public BaseMHDTimestepper<MHDModel, MessengerT>
{
    using Super = BaseMHDTimestepper<MHDModel, MessengerT>;

public:
    EulerIntegrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict, /*n_extra_states=*/0}
        , euler_{dict}
    {
    }

    void operator()(MHDModel& model, Super::MHDStateT& state,
                    Super::FluxT& fluxes, Super::Messenger& bc,
                    Super::level_t& level, double const currentTime,
                    double const newTime) override
    {
        this->resetButcherFluxes_(model, level);

        euler_(model, state, state, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state.E, fluxes, level);
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
    Euler<FVMethodStrategy, MHDModel> euler_;
};



} // namespace PHARE::solver

#endif
