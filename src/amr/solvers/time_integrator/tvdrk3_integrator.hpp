#ifndef PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_TVDRK3_INTEGRATOR_HPP

#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "amr/solvers/time_integrator/euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class TVDRK3Integrator
{
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
        : euler_{dict}
        , butcherFluxes_{{"timeRho_fx", core::MHDQuantity::Scalar::ScalarFlux_x},
                         {"timeRhoV_fx", core::MHDQuantity::Vector::VecFlux_x},
                         {"timeB_fx", core::MHDQuantity::Vector::VecFlux_x},
                         {"timeEtot_fx", core::MHDQuantity::Scalar::ScalarFlux_x},

                         {"timeRho_fy", core::MHDQuantity::Scalar::ScalarFlux_y},
                         {"timeRhoV_fy", core::MHDQuantity::Vector::VecFlux_y},
                         {"timeB_fy", core::MHDQuantity::Vector::VecFlux_y},
                         {"timeEtot_fy", core::MHDQuantity::Scalar::ScalarFlux_y},

                         {"timeRho_fz", core::MHDQuantity::Scalar::ScalarFlux_z},
                         {"timeRhoV_fz", core::MHDQuantity::Vector::VecFlux_z},
                         {"timeB_fz", core::MHDQuantity::Vector::VecFlux_z},
                         {"timeEtot_fz", core::MHDQuantity::Scalar::ScalarFlux_z}}
        , butcherE_{"timeE", core::MHDQuantity::Vector::E}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, level_t& level,
                    double const currentTime, double const newTime)
    {
        resetButcherFluxes_(model, level);

        // U1 = Euler(Un)
        euler_(model, state, state1_, fluxes, bc, level, currentTime, newTime);

        accumulateButcherFluxes_(model, state.E, fluxes, level, w01_ * w11_);

        // U1 = Euler(U1)
        euler_(model, state1_, state1_, fluxes, bc, level, currentTime, newTime);

        accumulateButcherFluxes_(model, state1_.E, fluxes, level, w01_ * w11_);

        // U2 = 0.75*Un + 0.25*U1
        tvdrk3_step_(level, model, newTime, state2_, RKPair_t{w00_, state},
                     RKPair_t{w01_, state1_});

        // U2 = Euler(U2)
        euler_(model, state2_, state2_, fluxes, bc, level, currentTime, newTime);

        accumulateButcherFluxes_(model, state2_.E, fluxes, level, w11_);

        // Un+1 = 1/3*Un + 2/3*Euler(U2)
        tvdrk3_step_(level, model, newTime, state, RKPair_t{w10_, state}, RKPair_t{w11_, state2_});
    }

    void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(state1_);
        model.resourcesManager->registerResources(state2_);
        model.resourcesManager->registerResources(butcherFluxes_);
        model.resourcesManager->registerResources(butcherE_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(state1_, patch, allocateTime);
        model.resourcesManager->allocate(state2_, patch, allocateTime);
        model.resourcesManager->allocate(butcherFluxes_, patch, allocateTime);
        model.resourcesManager->allocate(butcherE_, patch, allocateTime);
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
        return std::forward_as_tuple(state1_, state2_, butcherFluxes_, butcherE_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state1_, state2_, butcherFluxes_, butcherE_);
    }

    auto exposeFluxes() { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

    auto const exposeFluxes() const { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

private:
    void resetButcherFluxes_(MHDModel& model, auto& level)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _ = model.resourcesManager->setOnPatch(*patch, butcherFluxes_, butcherE_);

            evalFluxesOnGhostBox(
                layout, [&](auto& left, auto const&... args) mutable { left(args...) = 0.0; },
                butcherFluxes_);

            layout.evalOnGhostBox(butcherE_(core::Component::X), [&](auto const&... args) mutable {
                butcherE_(core::Component::X)(args...) = 0.0;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Y), [&](auto const&... args) mutable {
                butcherE_(core::Component::Y)(args...) = 0.0;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Z), [&](auto const&... args) mutable {
                butcherE_(core::Component::Z)(args...) = 0.0;
            });
        }
    }

    void accumulateButcherFluxes_(MHDModel& model, VecFieldT& E, auto& fluxes, auto& level,
                                  double const coef)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _
                = model.resourcesManager->setOnPatch(*patch, butcherFluxes_, butcherE_, fluxes, E);

            evalFluxesOnGhostBox(
                layout,
                [&](auto& left, auto const& right, auto const&... args) mutable {
                    if (std::isnan(left(args...)))
                    {
                        left(args...) = 0.0;
                    }
                    left(args...) += right(args...) * coef;
                },
                butcherFluxes_, fluxes);

            layout.evalOnGhostBox(butcherE_(core::Component::X), [&](auto const&... args) mutable {
                if (std::isnan(butcherE_(core::Component::X)(args...)))
                {
                    butcherE_(core::Component::X)(args...) = 0.0;
                }
                butcherE_(core::Component::X)(args...) += E(core::Component::X)(args...) * coef;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Y), [&](auto const&... args) mutable {
                if (std::isnan(butcherE_(core::Component::Y)(args...)))
                {
                    butcherE_(core::Component::Y)(args...) = 0.0;
                }
                butcherE_(core::Component::Y)(args...) += E(core::Component::Y)(args...) * coef;
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Z), [&](auto const&... args) mutable {
                if (std::isnan(butcherE_(core::Component::Z)(args...)))
                {
                    butcherE_(core::Component::Z)(args...) = 0.0;
                }
                butcherE_(core::Component::Z)(args...) += E(core::Component::Z)(args...) * coef;
            });
        }
    }

    static constexpr auto w00_{0.75};
    static constexpr auto w01_{0.25};
    static constexpr auto w10_{1. / 3.};
    static constexpr auto w11_{2. / 3.};

    Euler<FVMethodStrategy, MHDModel> euler_;
    core::AllFluxes<FieldT, VecFieldT> butcherFluxes_;
    VecFieldT butcherE_;
    RKUtils_t tvdrk3_step_;

    MHDStateT state1_{"state1"};
    MHDStateT state2_{"state2"};
};

} // namespace PHARE::solver

#endif
