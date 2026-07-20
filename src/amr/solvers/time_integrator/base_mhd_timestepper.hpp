#ifndef PHARE_CORE_NUMERICS_BASE_MHD_TIMESTEPPER_HPP
#define PHARE_CORE_NUMERICS_BASE_MHD_TIMESTEPPER_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/utilities/algorithm.hpp"

#include "initializer/data_provider.hpp"

#include <SAMRAI/hier/Patch.h>
#include <vector>

namespace PHARE::amr
{
template<typename MHDModel>
class MHDMessenger;
} // namespace PHARE::amr

namespace PHARE::solver
{
template<typename MHDModel, typename MessengerT = amr::MHDMessenger<MHDModel>>
class BaseMHDTimestepper
{
    using FieldT      = MHDModel::field_type;
    using VecFieldT   = MHDModel::vecfield_type;
    using GridLayoutT = MHDModel::gridlayout_type;

public:
    using MHDStateT = MHDModel::state_type;
    using FluxT     = core::AllFluxes<FieldT, VecFieldT>;
    using level_t   = MHDModel::level_t;
    using Messenger = MessengerT;

    BaseMHDTimestepper(PHARE::initializer::PHAREDict const& dict, std::size_t const n_extra_states)
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
        for (std::size_t i = 0; i < n_extra_states; ++i)
            extra_states_.emplace_back("state" + std::to_string(i + 1));
    }

    virtual ~BaseMHDTimestepper() = default;

    virtual void operator()(MHDModel& model, MHDStateT& state, FluxT& fluxes, Messenger& bc,
                            level_t& level, double currentTime, double newTime)
        = 0;

    virtual void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(*this);
    }

    virtual void allocate(MHDModel& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const
    {
        model.resourcesManager->allocate(*this, patch, allocateTime);
    }

    void fillMessengerInfo(auto& info) const
    {
        for (auto const& s : extra_states_)
        {
            info.ghostDensity.push_back(s.rho.name());
            info.ghostMomentum.push_back(s.rhoV.name());
            info.ghostTotalEnergy.push_back(s.Etot.name());
            info.ghostElectric.push_back(s.E.name());
            info.ghostMagnetic.push_back(s.B.name());
            info.ghostCurrent.push_back(s.J.name());
        }
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(butcherFluxes_, butcherE_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(butcherFluxes_, butcherE_);
    }

    NO_DISCARD std::vector<MHDStateT>& getRunTimeResourcesViewList() { return extra_states_; }

    NO_DISCARD std::vector<MHDStateT> const& getRunTimeResourcesViewList() const
    {
        return extra_states_;
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

    std::vector<MHDStateT> extra_states_;
    core::AllFluxes<FieldT, VecFieldT> butcherFluxes_;
    VecFieldT butcherE_;
};


} // namespace PHARE::solver

#endif
