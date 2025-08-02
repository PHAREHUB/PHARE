#ifndef PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP

#include "core/numerics/constrained_transport/constrained_transport.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/euler.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class EulerIntegrator
{
    using FieldT      = typename MHDModel::field_type;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayoutT = typename MHDModel::gridlayout_type;

public:
    EulerIntegrator(PHARE::initializer::PHAREDict const& dict)
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
    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, auto& level,
                    double const currentTime, double const newTime)
    {
        resetButcherFluxes_(model, level);

        euler_(model, state, state, fluxes, bc, level, currentTime, newTime);

        accumulateButcherFluxes_(model, fluxes, level);
    }

    void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(butcherFluxes_);
        model.resourcesManager->registerResources(butcherE_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(butcherFluxes_, patch, allocateTime);
        model.resourcesManager->allocate(butcherE_, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const {}

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(butcherFluxes_, butcherE_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(butcherFluxes_, butcherE_);
    }

    auto exposeFluxes() { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

    auto exposeFluxes() const { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

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

    void accumulateButcherFluxes_(MHDModel& model, auto& fluxes, auto& level)
    {
        for (auto& patch : level)
        {
            auto const& layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _ = model.resourcesManager->setOnPatch(*patch, butcherFluxes_, butcherE_, fluxes,
                                                        model.state.E);

            evalFluxesOnGhostBox(
                layout,
                [&](auto& left, auto const& right, auto const&... args) mutable {
                    left(args...) += right(args...);
                },
                butcherFluxes_, fluxes);


            layout.evalOnGhostBox(butcherE_(core::Component::X), [&](auto const&... args) mutable {
                butcherE_(core::Component::X)(args...)
                    += model.state.E(core::Component::X)(args...);
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Y), [&](auto const&... args) mutable {
                butcherE_(core::Component::Y)(args...)
                    += model.state.E(core::Component::Y)(args...);
            });

            layout.evalOnGhostBox(butcherE_(core::Component::Z), [&](auto const&... args) mutable {
                butcherE_(core::Component::Z)(args...)
                    += model.state.E(core::Component::Z)(args...);
            });
        }
    }

    Euler<FVMethodStrategy, MHDModel> euler_;
    core::AllFluxes<FieldT, VecFieldT> butcherFluxes_;
    VecFieldT butcherE_;
};
} // namespace PHARE::solver

#endif
