#ifndef PHARE_CORE_NUMERICS_BASE_MHD_TIMESTEPPER_HPP
#define PHARE_CORE_NUMERICS_BASE_MHD_TIMESTEPPER_HPP

#include "core/utilities/algorithm.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"

#include "initializer/data_provider.hpp"

#include "amr/resources_manager/amr_utils.hpp"

namespace PHARE::solver
{
template<typename MHDModel>
class BaseMHDTimestepper
{
    using FieldT      = typename MHDModel::field_type;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayoutT = typename MHDModel::gridlayout_type;

public:
    BaseMHDTimestepper(PHARE::initializer::PHAREDict const& dict)
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

protected:
    void resetButcherFluxes_(MHDModel& model, auto& level)
    {
        auto& rm = *model.resourcesManager;
        for (auto& _ : rm.enumerate(level, butcherFluxes_, butcherE_))
        {
            butcherFluxes_.zero();
            butcherE_.zero();
        }
    }

    void accumulateButcherFluxes_(MHDModel& model, auto& E, auto& fluxes, auto& level,
                                  double const coef = 1.0)
    {
        auto& rm = *model.resourcesManager;
        for (auto& _ : rm.enumerate(level, butcherFluxes_, butcherE_, fluxes, E))
        {
            core::operate<core::PlusEqualsProduct>(butcherFluxes_, fluxes, coef);
            core::operate<core::PlusEqualsProduct>(butcherE_, E, coef);
        }
    }

    core::AllFluxes<FieldT, VecFieldT> butcherFluxes_;
    VecFieldT butcherE_;
};


} // namespace PHARE::solver

#endif
