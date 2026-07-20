#ifndef PHARE_CORE_NUMERICS_SSPRK4_5_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_SSPRK4_5_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"
#include "amr/solvers/time_integrator/compute_fluxes.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"
#include "amr/solvers/solver_mhd_field_evolvers.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<typename FVMethodStrategy, typename MHDModel,
         typename MessengerT = amr::MHDMessenger<MHDModel>>
class SSPRK4_5Integrator : public BaseMHDTimestepper<MHDModel, MessengerT>
{
    using Super = BaseMHDTimestepper<MHDModel, MessengerT>;

    using VecFieldT = MHDModel::vecfield_type;
    using MHDStateT = MHDModel::state_type;

    using Dispatchers_t = Dispatchers<MHDModel>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

public:
    SSPRK4_5Integrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict, /*n_extra_states=*/4}
        , euler_{dict}
        , compute_fluxes_{dict}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, Super::MHDStateT& state, Super::FluxT& fluxes,
                    Super::Messenger& bc, Super::level_t& level, double const currentTime,
                    double const newTime) override
    {
        auto& state1 = this->extra_states_[0];
        auto& state2 = this->extra_states_[1];
        auto& state3 = this->extra_states_[2];
        auto& state4 = this->extra_states_[3];

        this->resetButcherFluxes_(model, level);

        auto const dt = newTime - currentTime;

        // U1 = Un + w0_*dt*F(Un)
        euler_(model, state, state1, fluxes, bc, level, currentTime, newTime, w0_ * dt);

        this->accumulateButcherFluxes_(
            model, state.E, fluxes, level,
            (w0_ * w11_ * w21_ * w31_ * w43_ + w0_ * w11_ * w21_ * w41_ + w0_ * w11_ * w40_));

        // U2 = w10_*Un + w11_*U1 + w12_*dt*F(U1)
        //
        // U2 = w10_*Un + w11_*U1
        RKUtils_t{level, model}(newTime, state2, RKPair_t{w10_, state}, RKPair_t{w11_, state1});

        // U2 = U2 + w12_*dt*F(U1)
        compute_fluxes_(model, state1, fluxes, bc, level, newTime);

        euler_using_butcher_fluxes_(model, state2, state2, state1.E, fluxes, bc, level, newTime,
                                    w12_ * dt);

        this->accumulateButcherFluxes_(
            model, state1.E, fluxes, level,
            (w12_ * w21_ * w31_ * w43_ + w12_ * w21_ * w41_ + w12_ * w40_));

        // U3 = w20_*Un + w21_*U2 + w22_*dt*F(U2)
        //
        // U3 = w20_*Un + w21_*U2
        RKUtils_t{level, model}(newTime, state3, RKPair_t{w20_, state}, RKPair_t{w21_, state2});

        // U3 = U3 + w22_*dt*F(U2)
        compute_fluxes_(model, state2, fluxes, bc, level, newTime);

        euler_using_butcher_fluxes_(model, state3, state3, state2.E, fluxes, bc, level, newTime,
                                    w22_ * dt);

        this->accumulateButcherFluxes_(model, state2.E, fluxes, level,
                                       (w22_ * w31_ * w43_ + w22_ * w41_));

        // U4 = w30_*Un + w31_*U3 + w32_*dt*F(U3)
        //
        // U4 = w30_*Un + w31_*U3
        RKUtils_t{level, model}(newTime, state4, RKPair_t{w30_, state}, RKPair_t{w31_, state3});

        // U4 = U4 + w32_*dt*F(U3)
        // if we were not using butcher formulation, we would need a separate flux buffer for F(U3)
        // for the final step
        compute_fluxes_(model, state3, fluxes, bc, level, newTime);

        euler_using_butcher_fluxes_(model, state4, state4, state3.E, fluxes, bc, level, newTime,
                                    w32_ * dt);

        this->accumulateButcherFluxes_(model, state3.E, fluxes, level, (w32_ * w43_ + w42_));

        compute_fluxes_(model, state4, fluxes, bc, level, newTime);

        this->accumulateButcherFluxes_(model, state4.E, fluxes, level, w44_);

        euler_using_butcher_fluxes_(model, state, state, this->butcherE_, this->butcherFluxes_, bc,
                                    level, newTime, dt);
    }

    void registerResources(MHDModel& model) override
    {
        Super::registerResources(model);
        euler_.registerResources(model);
        // no need to register compute_fluxes_, as it shares its resources with euler_
    }

    void allocate(MHDModel& model, SAMRAI::hier::Patch& patch,
                  double const allocateTime) const override
    {
        Super::allocate(model, patch, allocateTime);
        euler_.allocate(model, patch, allocateTime);
        // no need to allocate compute_fluxes_, as it shares its resources with euler_
    }

    using Super::exposeFluxes;

private:
    static constexpr auto w0_{0.391752226571890};
    static constexpr auto w10_{0.444370493651235};
    static constexpr auto w11_{0.555629506348765};
    static constexpr auto w12_{0.368410593050371};
    static constexpr auto w20_{0.620101851488403};
    static constexpr auto w21_{0.379898148511597};
    static constexpr auto w22_{0.251891774271694};
    static constexpr auto w30_{0.178079954393132};
    static constexpr auto w31_{0.821920045606868};
    static constexpr auto w32_{0.544974750228521};
    static constexpr auto w40_{0.517231671970585};
    static constexpr auto w41_{0.096059710526147};
    static constexpr auto w42_{0.063692468666290};
    static constexpr auto w43_{0.386708617503268};
    static constexpr auto w44_{0.226007483236906};

    Euler<FVMethodStrategy, MHDModel> euler_;
    ComputeFluxes<FVMethodStrategy, MHDModel> compute_fluxes_;
    EulerUsingComputedFlux<MHDModel> euler_using_butcher_fluxes_;
};

} // namespace PHARE::solver

#endif
