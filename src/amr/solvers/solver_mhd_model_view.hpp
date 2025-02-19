#ifndef PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP
#define PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP

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
    template<typename Layouts, typename States>
    void operator()(Layouts const& layouts, States& states)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], to_conservative_);
            to_conservative_(states.rho, states.V, states.B, states.P, states.rhoV, states.Etot);
        }
    }

    core_type to_conservative_;
};

template<typename GridLayout>
class ToPrimitiveTransformer
{
    using core_type = PHARE::core::ToPrimitiveConverter<GridLayout>;

public:
    template<typename Layouts, typename States>
    void operator()(Layouts const& layouts, States& states)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], to_primtitve_);
            to_primtitve_(states.rho, states.rhoV, states.B, states.Etot, states.V, states.P);
        }
    }

    core_type to_primtitve_;
};

template<typename GridLayout, template<typename> typename FVMethod>
class FVMethodTransformer
{
    using core_type = FVMethod<GridLayout>;

public:
    template<typename Layouts, typename StateViews, typename Fluxes>
    void operator()(Layouts const& layouts, StateViews& states, Fluxes& fluxes)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], fvm_);
            fvm_(states, fluxes);
        }
    }

    core_type fvm_;
};


template<typename GridLayout>
class FiniteVolumeEulerTransformer
{
    using core_type = PHARE::core::FiniteVolumeEuler<GridLayout>;

public:
    template<typename Layouts, typename StateViews, typename Fluxes>
    void operator()(Layouts const& layouts, StateViews const& states, StateViews statesnew,
                    double const dt, Fluxes const& fluxes)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], euler_);
            euler_(states, statesnew, dt, fluxes);
        }
    }

    core_type euler_;
};

template<typename GridLayout>
class ConstrainedTransportTransformer
{
    using core_type = PHARE::core::ConstrainedTransport<GridLayout>;

public:
    template<typename Layout, typename StateViews, typename Fluxes>
    void operator()(Layout const& layouts, StateViews& states, Fluxes const& fluxes)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], constrained_transport_);
            constrained_transport_(states.E, fluxes);
        }
    }

    core_type constrained_transport_;
};

template<typename GridLayout>
class FaradayMHDTransformer
{
    using core_type = PHARE::core::Faraday<GridLayout>;

public:
    template<typename GridLayouts, typename StateViews>
    void operator()(GridLayouts const& layouts, StateViews const& states, StateViews& statesnew,
                    double dt)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], faraday_);
            faraday_(states.B, states.E, statesnew.B, dt);
        }
    }

    core_type faraday_;
};

template<typename GridLayout>
class RKUtilsTransformer
{
    using core_type = PHARE::core::RKUtils<GridLayout>;

public:
    template<typename Layouts, typename ReturnState, typename... Pairs>
    void operator()(Layouts const& layouts, ReturnState& res, Pairs... pairs)
    {
        for (std::size_t i = 0; i < layouts.size(); ++i)
        {
            auto _ = core::SetLayout(layouts[i], rkutils_);
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

}; // namespace PHARE::solver

#endif // PHARE_SOLVER_SOLVER_MHD_MODEL_VIEW_HPP
