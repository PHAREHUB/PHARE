#ifndef PHARE_MHD_DEFAULTS_HPP
#define PHARE_MHD_DEFAULTS_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE
{
template<template<typename> typename FVmethod, typename MHDModel>
struct DefaultTimeIntegrator
{
    DefaultTimeIntegrator(PHARE::initializer::PHAREDict const& /*dict*/)
        : butcherFluxes_{{"timeRho_fx", core::MHDQuantity::Scalar::ScalarFlux_x},
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

    void operator()(MHDModel& /*model*/, MHDModel::state_type& /*state*/, auto& /*fluxes*/,
                    auto& /*fromCoarser*/, auto& /*level*/, double const /*currentTime*/,
                    double const /*newTime*/)
    {
    }

    void registerResources(MHDModel& /*model*/) {}

    void allocate(MHDModel& /*model*/, auto& /*patch*/, double const /*allocateTime*/) const {}

    void fillMessengerInfo(auto& /*info*/) const {}

    auto exposeFluxes() { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

    auto exposeFluxes() const { return std::forward_as_tuple(butcherFluxes_, butcherE_); }

    core::AllFluxes<typename MHDModel::field_type, typename MHDModel::vecfield_type> butcherFluxes_;
    MHDModel::vecfield_type butcherE_;
};

template<typename GridLayout, typename SlopeLimiter>
struct DefaultReconstruction
{
};

template<typename GridLayout, bool HallFlag>
struct DefaultRiemannSolver
{
};

template<bool Hall, bool Resistivity, bool HyperResistivity>
struct DefaultEquations
{
};
} // namespace PHARE

#endif
