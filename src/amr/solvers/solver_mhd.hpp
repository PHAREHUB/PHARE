#ifndef PHARE_SOLVER_MHD_HPP
#define PHARE_SOLVER_MHD_HPP

#include <array>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/mhd_messenger.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"

namespace PHARE::solver
{
template<typename MHDModel, typename AMR_Types, typename Messenger = amr::MHDMessenger<MHDModel>,
         typename ModelViews_t = MHDModelView<MHDModel>>
class SolverMHD : public ISolver<AMR_Types>
{
private:
    static constexpr auto dimension = MHDModel::dimension;

    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

    using FieldT      = typename MHDModel::field_type;
    using VecFieldT   = typename MHDModel::vecfield_type;
    using GridLayout  = typename MHDModel::gridlayout_type;
    using MHDQuantity = core::MHDQuantity;

    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using Direction        = core::Direction;

    using Ampere_t                  = typename MHDModelView<MHDModel>::Ampere_t;
    using GodunovFluxes_t           = typename MHDModelView<MHDModel>::GodunovFluxes_t;
    using ToPrimitiveConverter_t    = typename MHDModelView<MHDModel>::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = typename MHDModelView<MHDModel>::ToConservativeConverter_t;
    using TimeIntegrator_t          = typename MHDModelView<MHDModel>::TimeIntegrator_t;


    // Flux calculations
    FieldT rho_x{"rho_x", MHDQuantity::Scalar::ScalarFlux_x};
    VecFieldT rhoV_x{"rhoV_x", MHDQuantity::Vector::VecFlux_x};
    VecFieldT B_x{"B_x", MHDQuantity::Vector::VecFlux_x};
    FieldT Etot_x{"rho_x", MHDQuantity::Scalar::ScalarFlux_x};

    FieldT rho_y{"rho_y", MHDQuantity::Scalar::ScalarFlux_y};
    VecFieldT rhoV_y{"rhoV_y", MHDQuantity::Vector::VecFlux_y};
    VecFieldT B_y{"B_y", MHDQuantity::Vector::VecFlux_y};
    FieldT Etot_y{"rho_y", MHDQuantity::Scalar::ScalarFlux_y};

    FieldT rho_z{"rho_z", MHDQuantity::Scalar::ScalarFlux_z};
    VecFieldT rhoV_z{"rhoV_z", MHDQuantity::Vector::VecFlux_z};
    VecFieldT B_z{"B_z", MHDQuantity::Vector::VecFlux_z};
    FieldT Etot_z{"rho_z", MHDQuantity::Scalar::ScalarFlux_z};

    // Time integration
    FieldT rho1{"rho1", MHDQuantity::Scalar::rho};
    VecFieldT rhoV1{"rhoV1", MHDQuantity::Vector::rhoV};
    VecFieldT B1{"B1", MHDQuantity::Vector::B};
    FieldT Etot1{"Etot1", MHDQuantity::Scalar::Etot};

    FieldT rho2{"rho2", MHDQuantity::Scalar::rho};
    VecFieldT rhoV2{"rhoV2", MHDQuantity::Vector::rhoV};
    VecFieldT B2{"B2", MHDQuantity::Vector::B};
    FieldT Etot2{"Etot2", MHDQuantity::Scalar::Etot};

    Ampere_t ampere_;
    GodunovFluxes_t godunov_;
    ToPrimitiveConverter_t to_primitive_;
    ToConservativeConverter_t to_conservative_;
    TimeIntegrator_t time_integrator_;

    std::string const integrator_;

public:
    SolverMHD(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"MHDSolver"}
        , godunov_{dict["godunov"]}
        , to_primitive_{dict["to_primitive"]}
        , to_conservative_{dict["to_conservative"]}
        , integrator_{dict["integrator"].template to<std::string>()}
    {
    }

    virtual ~SolverMHD() = default;

    std::string modelName() const override { return MHDModel::model_name; }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override {}

    void registerResources(IPhysicalModel<AMR_Types>& /*model*/) override {}

    // TODO make this a resourcesUser
    void allocate(IPhysicalModel<AMR_Types>& /*model*/, patch_t& /*patch*/,
                  double const /*allocateTime*/) const override
    {
    }

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& view,
                      IMessenger& fromCoarserMessenger, const double currentTime,
                      const double newTime) override;


    std::shared_ptr<ISolverModelView> make_view(level_t& level, IPhysicalModel_t& model) override
    {
        /*return std::make_shared<ModelViews_t>(level, dynamic_cast<MHDModel&>(model));*/
        throw std::runtime_error("no MHD model yet");
        return nullptr;
    }

private:
    void godunov_fluxes_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                         double const currentTime, double const newTime);

    void time_integration_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                           double const currentTime, double const newTime);

    template<typename Layout, typename View, typename VecField, typename Qs, typename... Fluxes>
    void euler_(Messenger& fromCoarser, level_t& level, double const newTime, Layout& layouts,
                Qs quantities, View E, VecField Emodel, Qs quantities_new, double const dt,
                Fluxes... fluxes);

    struct TimeSetter
    {
        /*template <typename QuantityAccessor>*/
        /*void operator()(QuantityAccessor accessor) {*/
        /*    for (auto& state : views)*/
        /*        views.model().resourcesManager->setTime(accessor(state), *state.patch, newTime);*/
        /*}*/
        /**/
        /*ModelViews_t& views;*/
        /*double        newTime;*/
    };
};

// -----------------------------------------------------------------------------

template<typename MHDModel, typename AMR_Types, typename Messenger, typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, Messenger, ModelViews_t>::advanceLevel(
    hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& view,
    IMessenger& fromCoarserMessenger, const double currentTime, const double newTime)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::advanceLevel");

    auto& modelView   = dynamic_cast<ModelViews_t&>(view);
    auto& fromCoarser = dynamic_cast<Messenger&>(fromCoarserMessenger);
    auto level        = hierarchy.getPatchLevel(levelNumber);

    godunov_fluxes_(*level, modelView, fromCoarser, currentTime, newTime);

    time_integration_(*level, modelView, fromCoarser, currentTime, newTime);
}


template<typename MHDModel, typename AMR_Types, typename Messenger, typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, Messenger, ModelViews_t>::godunov_fluxes_(
    level_t& level, ModelViews_t& views, Messenger& fromCoarser, double const currentTime,
    double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::godunov_fluxes_");

    fromCoarser.fillMagneticGhosts(views.model().state.B, level, newTime);

    ampere_(views.layouts, views.B, views.J);

    fromCoarser.fillMomentGhosts(views.model().state.rho, level, newTime);
    fromCoarser.fillMomentGhosts(views.model().state.V(core::Component::X), level, newTime);
    fromCoarser.fillMomentGhosts(views.model().state.V(core::Component::Y), level, newTime);
    fromCoarser.fillMomentGhosts(views.model().state.V(core::Component::Z), level, newTime);
    fromCoarser.fillMomentGhosts(views.model().state.P, level, newTime);

    fromCoarser.fillCurrentGhosts(views.model().state.J, level, newTime);

    if constexpr (dimension == 1)
    {
        godunov_(views.layouts, views.rho, views.V, views.B, views.P, views.J, views.rho_x,
                 views.rhoV_x, views.B_x, views.Etot_x);
    }
    if constexpr (dimension == 2)
    {
        godunov_(views.layouts, views.rho, views.V, views.B, views.P, views.J, views.rho_x,
                 views.rhoV_x, views.B_x, views.Etot_x, views.rho_y, views.rhoV_y, views.B_y,
                 views.Etot_y);
    }
    if constexpr (dimension == 3)
    {
        godunov_(views.layouts, views.rho, views.V, views.B, views.P, views.J, views.rho_x,
                 views.rhoV_x, views.B_x, views.Etot_x, views.rho_y, views.rhoV_y, views.B_y,
                 views.Etot_y, views.rho_z, views.rhoV_z, views.B_z, views.Etot_z);
    }
}


template<typename MHDModel, typename AMR_Types, typename Messenger, typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, Messenger, ModelViews_t>::time_integration_(
    level_t& level, ModelViews_t& views, Messenger& fromCoarser, double const currentTime,
    double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::time_integrator_");

    auto dt = newTime - currentTime;

    to_conservative_(views.layouts, views.rho, views.V, views.B, views.P, views.rhoV, views.Etot);

    fromCoarser.fillMagneticFluxGhosts(views.B_x, level, newTime);

    if constexpr (dimension == 1)
    {
        if (integrator_ == "euler")
            time_integrator_.euler(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                   views.E, dt, views.rho_x, views.rhoV_x, views.B_x, views.Etot_x);
        else if (integrator_ == "tvdrk2")
            time_integrator_.tvdrk2(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                    views.rho1, views.rhoV1, views.B1, views.Etot1, views.E, dt,
                                    views.rho_x, views.rhoV_x, views.B_x, views.Etot_x);
        else if (integrator_ == "tvdrk3")
            time_integrator_.tvdrk3(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                    views.rho1, views.rhoV1, views.B1, views.Etot1, views.rho2,
                                    views.rhoV2, views.B2, views.Etot2, views.E, dt, views.rho_x,
                                    views.rhoV_x, views.B_x, views.Etot_x);
    }
    if constexpr (dimension >= 2)
    {
        fromCoarser.fillMagneticFluxGhosts(views.B_y, level, newTime);

        if constexpr (dimension == 2)
        {
            if (integrator_ == "euler")
                time_integrator_.euler(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                       views.E, dt, views.rho_x, views.rhoV_x, views.B_x,
                                       views.Etot_x, views.rho_y, views.rhoV_y, views.B_y,
                                       views.Etot_y);
            else if (integrator_ == "tvdrk2")
                time_integrator_.tvdrk2(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                        views.rho1, views.rhoV1, views.B1, views.Etot1, views.E, dt,
                                        views.rho_x, views.rhoV_x, views.B_x, views.Etot_x,
                                        views.rho_y, views.rhoV_y, views.B_y, views.Etot_y);
            else if (integrator_ == "tvdrk3")
                time_integrator_.tvdrk3(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                        views.rho1, views.rhoV1, views.B1, views.Etot1, views.rho2,
                                        views.rhoV2, views.B2, views.Etot2, views.E, dt,
                                        views.rho_x, views.rhoV_x, views.B_x, views.Etot_x,
                                        views.rho_y, views.rhoV_y, views.B_y, views.Etot_y);
        }
        if constexpr (dimension == 3)
        {
            fromCoarser.fillMagneticFluxGhosts(views.B_z, level, newTime);

            if (integrator_ == "euler")
                time_integrator_.euler(
                    views.layouts, views.rho, views.rhoV, views.B, views.Etot, views.E, dt,
                    views.rho_x, views.rhoV_x, views.B_x, views.Etot_x, views.rho_y, views.rhoV_y,
                    views.B_y, views.Etot_y, views.rho_z, views.rhoV_z, views.B_z, views.Etot_z);
            else if (integrator_ == "tvdrk2")
                time_integrator_.tvdrk2(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                        views.rho1, views.rhoV1, views.B1, views.Etot1, views.E, dt,
                                        views.rho_x, views.rhoV_x, views.B_x, views.Etot_x,
                                        views.rho_y, views.rhoV_y, views.B_y, views.Etot_y,
                                        views.rho_z, views.rhoV_z, views.B_z, views.Etot_z);
            else if (integrator_ == "tvdrk3")
                time_integrator_.tvdrk3(views.layouts, views.rho, views.rhoV, views.B, views.Etot,
                                        views.rho1, views.rhoV1, views.B1, views.Etot1, views.rho2,
                                        views.rhoV2, views.B2, views.Etot2, views.E, dt,
                                        views.rho_x, views.rhoV_x, views.B_x, views.Etot_x,
                                        views.rho_y, views.rhoV_y, views.B_y, views.Etot_y,
                                        views.rho_z, views.rhoV_z, views.B_z, views.Etot_z);
        }
    }

    to_primitive_(views.layouts, views.rho, views.rhoV, views.B, views.Etot, views.V, views.P);
}

} // namespace PHARE::solver

#endif
