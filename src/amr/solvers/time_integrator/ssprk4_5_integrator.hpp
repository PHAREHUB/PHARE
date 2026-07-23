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
template<typename FVMethodStrategy, typename MHDModel>
class SSPRK4_5Integrator : public BaseMHDTimestepper<MHDModel>
{
    using Super = BaseMHDTimestepper<MHDModel>;

    using level_t     = typename MHDModel::level_t;
    using FieldT      = typename MHDModel::field_type;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayoutT = typename MHDModel::gridlayout_type;
    using MHDStateT   = typename MHDModel::state_type;

    using Dispatchers_t = Dispatchers<MHDModel>;
    using RKUtils_t     = Dispatchers_t::RKUtils_t;

    using RKPair_t = core::RKPair<typename VecFieldT::value_type, MHDStateT>;

public:
    SSPRK4_5Integrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
        , euler_{dict}
        , compute_fluxes_{dict}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, level_t& level,
                    double const currentTime, double const newTime)
    {
        this->resetButcherFluxes_(model, level);

        auto const dt = newTime - currentTime;

        // U1 = Un + w0_*dt*F(Un). Fill U1's coarse-fine ghosts at its abscissa t_n + c1_*dt.
        euler_(model, state, state1_, fluxes, bc, level, currentTime, currentTime + c1_ * dt,
               w0_ * dt);

        this->accumulateButcherFluxes_(
            model, state.E, fluxes, level,
            (w0_ * w11_ * w21_ * w31_ * w43_ + w0_ * w11_ * w21_ * w41_ + w0_ * w11_ * w40_));

        // U2 = w10_*Un + w11_*U1 + w12_*dt*F(U1)
        //
        // U2 = w10_Un + w11_*U1
        RKUtils_t{level, model}(state2_, RKPair_t{w10_, state}, RKPair_t{w11_, state1_});

        // U2 = U2 + w12_*dt*F(U1)
        compute_fluxes_(model, state1_, fluxes, bc, level, newTime);

        euler_using_butcher_fluxes_(model, state2_, state2_, state1_.E, fluxes, bc, level,
                                    currentTime + c2_ * dt, w12_ * dt);

        this->accumulateButcherFluxes_(
            model, state1_.E, fluxes, level,
            (w12_ * w21_ * w31_ * w43_ + w12_ * w21_ * w41_ + w12_ * w40_));

        // U3 = w20_*Un + w21_*U2 + w22_*dt*F(U2)
        //
        // U3 = w20_*Un + w21_*U2
        RKUtils_t{level, model}(state3_, RKPair_t{w20_, state}, RKPair_t{w21_, state2_});

        // U3 = U3 + w22_*dt*F(U2)
        compute_fluxes_(model, state2_, fluxes, bc, level, newTime);

        euler_using_butcher_fluxes_(model, state3_, state3_, state2_.E, fluxes, bc, level,
                                    currentTime + c3_ * dt, w22_ * dt);

        this->accumulateButcherFluxes_(model, state2_.E, fluxes, level,
                                       (w22_ * w31_ * w43_ + w22_ * w41_));

        // U4 = w30_*Un + w31_*U3 + w32_*dt*F(U3)
        //
        // U4 = w30_*Un + w31_*U3
        RKUtils_t{level, model}(state4_, RKPair_t{w30_, state}, RKPair_t{w31_, state3_});

        // U4 = U4 + w32_*dt*F(U3)
        // if we were not using butcher formulation, we would need a separate flux buffer for F(U3)
        // for the final step
        compute_fluxes_(model, state3_, fluxes, bc, level, newTime);

        euler_using_butcher_fluxes_(model, state4_, state4_, state3_.E, fluxes, bc, level,
                                    currentTime + c4_ * dt, w32_ * dt);

        this->accumulateButcherFluxes_(model, state3_.E, fluxes, level, (w32_ * w43_ + w42_));

        compute_fluxes_(model, state4_, fluxes, bc, level, newTime);

        this->accumulateButcherFluxes_(model, state4_.E, fluxes, level, w44_);

        euler_using_butcher_fluxes_(model, state, state, this->butcherE_, this->butcherFluxes_, bc,
                                    level, newTime, dt);

        // Un+1 = w40_*U2 + w41_*U3 + w42_*F(U3) + w43_*U4 + w44_*dt*F(U4)
    }

    void registerResources(MHDModel& model)
    {
        Super::registerResources(model);
        model.resourcesManager->registerResources(state1_);
        model.resourcesManager->registerResources(state2_);
        model.resourcesManager->registerResources(state3_);
        model.resourcesManager->registerResources(state4_);
        euler_.registerResources(model);
        // probably we should have the same resources for euler and compute_fluxes
        // compute_fluxes_.registerResources(model);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        Super::allocate(model, patch, allocateTime);
        model.resourcesManager->allocate(state1_, patch, allocateTime);
        model.resourcesManager->allocate(state2_, patch, allocateTime);
        model.resourcesManager->allocate(state3_, patch, allocateTime);
        model.resourcesManager->allocate(state4_, patch, allocateTime);
        euler_.allocate(model, patch, allocateTime);
        // probably we should have the same resources for euler and compute_fluxes
        // compute_fluxes_.allocate(model, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const
    {
        auto fill_info = [&](auto& state) {
            info.ghostDensity.push_back(state.rho.name());
            info.ghostMomentum.push_back(state.rhoV.name());
            info.ghostTotalEnergy.push_back(state.Etot.name());
            info.ghostElectric.push_back(state.E.name());
            info.ghostMagnetic.push_back(state.B.name());
        };

        fill_info(state1_);
        fill_info(state2_);
        fill_info(state3_);
        fill_info(state4_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::tuple_cat(Super::getCompileTimeResourcesViewList(),
                              std::forward_as_tuple(state1_, state2_, state3_, state4_));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::tuple_cat(Super::getCompileTimeResourcesViewList(),
                              std::forward_as_tuple(state1_, state2_, state3_, state4_));
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

    // Stage abscissae c_i (row sums of the Shu-Osher tableau): the physical time
    // t_n + c_i*dt that each intermediate node approximates. Used to fill coarse-fine ghosts
    // by linear time interpolation at the node's own time. Match Spiteri-Ruuth SSPRK(5,4)
    // c = [0.39175, 0.58608, 0.47454, 0.93501].
    static constexpr double c1_{w0_};                // state1_
    static constexpr double c2_{w11_ * c1_ + w12_};  // state2_
    static constexpr double c3_{w21_ * c2_ + w22_};  // state3_
    static constexpr double c4_{w31_ * c3_ + w32_};  // state4_

    Euler<FVMethodStrategy, MHDModel> euler_;
    ComputeFluxes<FVMethodStrategy, MHDModel> compute_fluxes_;
    EulerUsingComputedFlux<MHDModel> euler_using_butcher_fluxes_;

    MHDStateT state1_{"state1"};
    MHDStateT state2_{"state2"};
    MHDStateT state3_{"state3"};
    MHDStateT state4_{"state4"};
};

} // namespace PHARE::solver

#endif
