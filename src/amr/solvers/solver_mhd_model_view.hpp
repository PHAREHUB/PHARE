#ifndef PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP
#define PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP

#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver.hpp"
#include "core/numerics/constrained_transport/constrained_transport.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/time_integrator_utils.hpp"

namespace PHARE::solver
{
template<typename GridLayout>
class ToConservativeTransformer
{
    using core_type = PHARE::core::ToConservativeConverter<GridLayout>;

public:
    template<typename MHDModel>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& state)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, state);
            auto _sl    = core::SetLayout(&layout, to_conservative_);
            to_conservative_(state.rho, state.V, state.B, state.P, state.rhoV, state.Etot);
        }
    }

    core_type to_conservative_;
};

template<typename GridLayout>
class ToPrimitiveTransformer
{
    using core_type = PHARE::core::ToPrimitiveConverter<GridLayout>;

public:
    template<typename MHDModel>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& state)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, state);
            auto _sl    = core::SetLayout(&layout, to_primitive_);
            to_primitive_(state.rho, state.rhoV, state.B, state.Etot, state.V, state.P);
        }
    }

    core_type to_primitive_;
};

template<typename GridLayout, template<typename> typename FVMethod>
class FVMethodTransformer
{
    using core_type = FVMethod<GridLayout>;

public:
    template<typename MHDModel>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& state,
                    auto& fluxes)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, state, fluxes);
            auto _sl    = core::SetLayout(&layout, fvm_);
            fvm_(state, fluxes);
        }
    }

    core_type fvm_;
};


template<typename GridLayout>
class FiniteVolumeEulerTransformer
{
    using core_type = PHARE::core::FiniteVolumeEuler<GridLayout>;

public:
    template<typename MHDModel>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& state,
                    MHDModel::state_type& statenew, auto& fluxes, double const dt)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, state, statenew, fluxes);
            auto _sl    = core::SetLayout(&layout, euler_);
            euler_(state, statenew, fluxes, dt);
        }
    }

    core_type euler_;
};

template<typename GridLayout>
class ConstrainedTransportTransformer
{
    using core_type = PHARE::core::ConstrainedTransport<GridLayout>;

public:
    template<typename MHDModel>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& state,
                    auto& fluxes)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, state, fluxes);
            auto _sl    = core::SetLayout(&layout, constrained_transport_);
            constrained_transport_(state.E, fluxes);
        }
    }

    core_type constrained_transport_;
};

template<typename GridLayout>
class FaradayMHDTransformer
{
    using core_type = PHARE::core::Faraday<GridLayout>;

public:
    template<typename MHDModel>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& state,
                    MHDModel::state_type& statenew, double dt)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, state, statenew);
            auto _sl    = core::SetLayout(&layout, faraday_);
            faraday_(state.B, state.E, statenew.B, dt);
        }
    }

    core_type faraday_;
};

template<typename GridLayout>
class RKUtilsTransformer
{
    using core_type = PHARE::core::RKUtils<GridLayout>;

public:
    template<typename MHDModel, typename... Pairs>
    void operator()(MHDModel::level_t const& level, MHDModel& model, MHDModel::state_type& res,
                    Pairs... pairs)
    {
        for (auto const& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _sp    = model.resourcesManager->setOnPatch(*patch, res, pairs.state...);
            auto _sl    = core::SetLayout(&layout, rkutils_);
            rkutils_(res, pairs...);
        }
    }

    core_type rkutils_;
};


template<typename GridLayout>
class Dispatchers
{
public:
    using ToPrimitiveConverter_t    = ToPrimitiveTransformer<GridLayout>;
    using ToConservativeConverter_t = ToConservativeTransformer<GridLayout>;

    template<template<typename> typename FVMethodStrategy>
    using FVMethod_t = FVMethodTransformer<GridLayout, FVMethodStrategy>;

    using FiniteVolumeEuler_t    = FiniteVolumeEulerTransformer<GridLayout>;
    using ConstrainedTransport_t = ConstrainedTransportTransformer<GridLayout>;
    using Faraday_t              = FaradayMHDTransformer<GridLayout>;
    using RKUtils_t              = RKUtilsTransformer<GridLayout>;
};

// for now keep identical interface as hybrid for simplicity
template<typename MHDModel_>
class MHDModelView : public ISolverModelView
{
public:
    using MHDModel_t       = MHDModel_;
    using level_t          = typename MHDModel_t::level_t;
    using IPhysicalModel_t = MHDModel_t::Interface;

    MHDModelView(level_t& level, IPhysicalModel_t& model)
        : model_{dynamic_cast<MHDModel_&>(model)}
    {
    }

    auto& model() { return model_; }
    auto& model() const { return model_; }

    MHDModel_t& model_;
};

}; // namespace PHARE::solver

#endif // PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP
