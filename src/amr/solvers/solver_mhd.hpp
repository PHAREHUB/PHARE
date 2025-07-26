#ifndef PHARE_SOLVER_MHD_HPP
#define PHARE_SOLVER_MHD_HPP

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "initializer/data_provider.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/mhd_messenger.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::solver
{
template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy,
         typename Messenger    = amr::MHDMessenger<MHDModel>,
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
    using MHDStateT   = typename MHDModel::state_type;
    using GridLayout  = typename MHDModel::gridlayout_type;
    using MHDQuantity = core::MHDQuantity;

    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;

    core::AllFluxes<FieldT, VecFieldT> fluxes_;

    TimeIntegratorStrategy evolve_;

    // Refluxing
    MHDStateT stateOld_{this->name() + "_stateOld"};

    core::AllFluxes<FieldT, VecFieldT> fluxSum_;
    VecFieldT fluxSumE_{this->name() + "_fluxSumE", MHDQuantity::Vector::E};

    std::unordered_map<std::size_t, double> oldTime_;

public:
    SolverMHD(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"MHDSolver"}
        , fluxes_{{"rho_fx", MHDQuantity::Scalar::ScalarFlux_x},
                  {"rhoV_fx", MHDQuantity::Vector::VecFlux_x},
                  {"B_fx", MHDQuantity::Vector::VecFlux_x},
                  {"Etot_fx", MHDQuantity::Scalar::ScalarFlux_x},

                  {"rho_fy", MHDQuantity::Scalar::ScalarFlux_y},
                  {"rhoV_fy", MHDQuantity::Vector::VecFlux_y},
                  {"B_fy", MHDQuantity::Vector::VecFlux_y},
                  {"Etot_fy", MHDQuantity::Scalar::ScalarFlux_y},

                  {"rho_fz", MHDQuantity::Scalar::ScalarFlux_z},
                  {"rhoV_fz", MHDQuantity::Vector::VecFlux_z},
                  {"B_fz", MHDQuantity::Vector::VecFlux_z},
                  {"Etot_fz", MHDQuantity::Scalar::ScalarFlux_z}}
        , evolve_{dict}
        , fluxSum_{{"sumRho_fx", MHDQuantity::Scalar::ScalarFlux_x},
                   {"sumRhoV_fx", MHDQuantity::Vector::VecFlux_x},
                   {"sumB_fx", MHDQuantity::Vector::VecFlux_x},
                   {"sumEtot_fx", MHDQuantity::Scalar::ScalarFlux_x},

                   {"sumRho_fy", MHDQuantity::Scalar::ScalarFlux_y},
                   {"sumRhoV_fy", MHDQuantity::Vector::VecFlux_y},
                   {"sumB_fy", MHDQuantity::Vector::VecFlux_y},
                   {"sumEtot_fy", MHDQuantity::Scalar::ScalarFlux_y},

                   {"sumRho_fz", MHDQuantity::Scalar::ScalarFlux_z},
                   {"sumRhoV_fz", MHDQuantity::Vector::VecFlux_z},
                   {"sumB_fz", MHDQuantity::Vector::VecFlux_z},
                   {"sumEtot_fz", MHDQuantity::Scalar::ScalarFlux_z}}
    {
    }

    virtual ~SolverMHD() = default;

    std::string modelName() const override { return MHDModel::model_name; }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;

    void registerResources(IPhysicalModel<AMR_Types>& model) override;

    // TODO make this a resourcesUser
    void allocate(IPhysicalModel<AMR_Types>& model, patch_t& patch,
                  double const allocateTime) const override;

    void prepareStep(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                     double const currentTime) override;

    void accumulateFluxSum(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                           double const coef) override;

    void resetFluxSum(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level) override;

    void reflux(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, IMessenger& messenger,
                double const time) override;

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& view,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;

    void onRegrid() override {}

    std::shared_ptr<ISolverModelView> make_view(level_t& level, IPhysicalModel_t& model) override
    {
        return std::make_shared<ModelViews_t>(level, dynamic_cast<MHDModel&>(model));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(fluxes_, fluxSum_, fluxSumE_, stateOld_, evolve_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(fluxes_, fluxSum_, fluxSumE_, stateOld_, evolve_);
    }

private:
    struct TimeSetter
    {
        template<typename QuantityAccessor>
        void operator()(QuantityAccessor accessor)
        {
            for (auto& state : views)
                views.model().resourcesManager->setTime(accessor(state), *state.patch, newTime);
        }

        ModelViews_t& views;
        double newTime;
    };
};

// -----------------------------------------------------------------------------

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger,
               ModelViews_t>::registerResources(IPhysicalModel_t& model)
{
    auto& mhdmodel = dynamic_cast<MHDModel&>(model);

    mhdmodel.resourcesManager->registerResources(fluxes_.rho_fx);
    mhdmodel.resourcesManager->registerResources(fluxes_.rhoV_fx);
    mhdmodel.resourcesManager->registerResources(fluxes_.B_fx);
    mhdmodel.resourcesManager->registerResources(fluxes_.Etot_fx);

    if constexpr (dimension >= 2)
    {
        mhdmodel.resourcesManager->registerResources(fluxes_.rho_fy);
        mhdmodel.resourcesManager->registerResources(fluxes_.rhoV_fy);
        mhdmodel.resourcesManager->registerResources(fluxes_.B_fy);
        mhdmodel.resourcesManager->registerResources(fluxes_.Etot_fy);

        if constexpr (dimension == 3)
        {
            mhdmodel.resourcesManager->registerResources(fluxes_.rho_fz);
            mhdmodel.resourcesManager->registerResources(fluxes_.rhoV_fz);
            mhdmodel.resourcesManager->registerResources(fluxes_.B_fz);
            mhdmodel.resourcesManager->registerResources(fluxes_.Etot_fz);
        }
    }

    mhdmodel.resourcesManager->registerResources(fluxSum_.rho_fx);
    mhdmodel.resourcesManager->registerResources(fluxSum_.rhoV_fx);
    mhdmodel.resourcesManager->registerResources(fluxSum_.B_fx);
    mhdmodel.resourcesManager->registerResources(fluxSum_.Etot_fx);

    if constexpr (dimension >= 2)
    {
        mhdmodel.resourcesManager->registerResources(fluxSum_.rho_fy);
        mhdmodel.resourcesManager->registerResources(fluxSum_.rhoV_fy);
        mhdmodel.resourcesManager->registerResources(fluxSum_.B_fy);
        mhdmodel.resourcesManager->registerResources(fluxSum_.Etot_fy);

        if constexpr (dimension == 3)
        {
            mhdmodel.resourcesManager->registerResources(fluxSum_.rho_fz);
            mhdmodel.resourcesManager->registerResources(fluxSum_.rhoV_fz);
            mhdmodel.resourcesManager->registerResources(fluxSum_.B_fz);
            mhdmodel.resourcesManager->registerResources(fluxSum_.Etot_fz);
        }
    }
    mhdmodel.resourcesManager->registerResources(fluxSumE_);

    mhdmodel.resourcesManager->registerResources(stateOld_);

    evolve_.registerResources(mhdmodel);
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger, ModelViews_t>::allocate(
    IPhysicalModel_t& model, patch_t& patch, double const allocateTime) const

{
    auto& mhdmodel = dynamic_cast<MHDModel&>(model);

    mhdmodel.resourcesManager->allocate(fluxes_.rho_fx, patch, allocateTime);
    mhdmodel.resourcesManager->allocate(fluxes_.rhoV_fx, patch, allocateTime);
    mhdmodel.resourcesManager->allocate(fluxes_.B_fx, patch, allocateTime);
    mhdmodel.resourcesManager->allocate(fluxes_.Etot_fx, patch, allocateTime);

    if constexpr (dimension >= 2)
    {
        mhdmodel.resourcesManager->allocate(fluxes_.rho_fy, patch, allocateTime);
        mhdmodel.resourcesManager->allocate(fluxes_.rhoV_fy, patch, allocateTime);
        mhdmodel.resourcesManager->allocate(fluxes_.B_fy, patch, allocateTime);
        mhdmodel.resourcesManager->allocate(fluxes_.Etot_fy, patch, allocateTime);

        if constexpr (dimension == 3)
        {
            mhdmodel.resourcesManager->allocate(fluxes_.rho_fz, patch, allocateTime);
            mhdmodel.resourcesManager->allocate(fluxes_.rhoV_fz, patch, allocateTime);
            mhdmodel.resourcesManager->allocate(fluxes_.B_fz, patch, allocateTime);
            mhdmodel.resourcesManager->allocate(fluxes_.Etot_fz, patch, allocateTime);
        }
    }

    mhdmodel.resourcesManager->allocate(fluxSum_.rho_fx, patch, allocateTime);
    mhdmodel.resourcesManager->allocate(fluxSum_.rhoV_fx, patch, allocateTime);
    mhdmodel.resourcesManager->allocate(fluxSum_.B_fx, patch, allocateTime);
    mhdmodel.resourcesManager->allocate(fluxSum_.Etot_fx, patch, allocateTime);

    if constexpr (dimension >= 2)
    {
        mhdmodel.resourcesManager->allocate(fluxSum_.rho_fy, patch, allocateTime);
        mhdmodel.resourcesManager->allocate(fluxSum_.rhoV_fy, patch, allocateTime);
        mhdmodel.resourcesManager->allocate(fluxSum_.B_fy, patch, allocateTime);
        mhdmodel.resourcesManager->allocate(fluxSum_.Etot_fy, patch, allocateTime);

        if constexpr (dimension == 3)
        {
            mhdmodel.resourcesManager->allocate(fluxSum_.rho_fz, patch, allocateTime);
            mhdmodel.resourcesManager->allocate(fluxSum_.rhoV_fz, patch, allocateTime);
            mhdmodel.resourcesManager->allocate(fluxSum_.B_fz, patch, allocateTime);
            mhdmodel.resourcesManager->allocate(fluxSum_.Etot_fz, patch, allocateTime);
        }
    }
    mhdmodel.resourcesManager->allocate(fluxSumE_, patch, allocateTime);

    mhdmodel.resourcesManager->allocate(stateOld_, patch, allocateTime);

    evolve_.allocate(mhdmodel, patch, allocateTime);
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger,
               ModelViews_t>::fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info)
    const

{
    auto& mhdInfo = dynamic_cast<amr::MHDMessengerInfo&>(*info);

    mhdInfo.ghostMagneticFluxesX.emplace_back(core::VecFieldNames{fluxes_.B_fx});

    if constexpr (dimension >= 2)
    {
        mhdInfo.ghostMagneticFluxesY.emplace_back(core::VecFieldNames{fluxes_.B_fy});

        if constexpr (dimension == 3)
        {
            mhdInfo.ghostMagneticFluxesZ.emplace_back(core::VecFieldNames{fluxes_.B_fz});
        }
    }

    evolve_.fillMessengerInfo(mhdInfo);

    auto&& [timeFluxes, timeElectric] = evolve_.exposeFluxes();

    mhdInfo.reflux          = core::AllFluxesNames{timeFluxes};
    mhdInfo.refluxElectric  = core::VecFieldNames{timeElectric};
    mhdInfo.fluxSum         = core::AllFluxesNames{fluxSum_};
    mhdInfo.fluxSumElectric = core::VecFieldNames{fluxSumE_};

    // for the faraday in reflux
    mhdInfo.ghostElectric.emplace_back(core::VecFieldNames{timeElectric});
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger, ModelViews_t>::prepareStep(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, double const currentTime)
{
    oldTime_[level.getLevelNumber()] = currentTime;

    auto& mhdModel = dynamic_cast<MHDModel&>(model);

    auto& rho  = mhdModel.state.rho;
    auto& rhoV = mhdModel.state.rhoV;
    auto& B    = mhdModel.state.B;
    auto& Etot = mhdModel.state.Etot;

    for (auto& patch : level)
    {
        auto dataOnPatch
            = mhdModel.resourcesManager->setOnPatch(*patch, rho, rhoV, B, Etot, stateOld_);

        mhdModel.resourcesManager->setTime(stateOld_.rho, *patch, currentTime);
        mhdModel.resourcesManager->setTime(stateOld_.rhoV, *patch, currentTime);
        mhdModel.resourcesManager->setTime(stateOld_.B, *patch, currentTime);
        mhdModel.resourcesManager->setTime(stateOld_.Etot, *patch, currentTime);

        stateOld_.rho.copyData(rho);
        stateOld_.rhoV.copyData(rhoV);
        stateOld_.B.copyData(B);
        stateOld_.Etot.copyData(Etot);
    }
}


template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger,
               ModelViews_t>::accumulateFluxSum(IPhysicalModel_t& model,
                                                SAMRAI::hier::PatchLevel& level, double const coef)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::accumulateFluxSum");

    auto& mhdModel = dynamic_cast<MHDModel&>(model);

    for (auto& patch : level)
    {
        auto&& [timeFluxes, timeElectric] = evolve_.exposeFluxes();
        auto const& layout                = amr::layoutFromPatch<GridLayout>(*patch);
        auto _ = mhdModel.resourcesManager->setOnPatch(*patch, fluxSum_, fluxSumE_, timeFluxes,
                                                       timeElectric);

        evalFluxesOnGhostBox(
            layout,
            [&](auto& left, auto const& right, auto const&... args) mutable {
                if (std::isnan(left(args...)))
                {
                    left(args...) = 0.0;
                }
                left(args...) += right(args...) * coef; // fluxSum_ += exposedFluxes*coef
            },
            fluxSum_, timeFluxes);

        layout.evalOnGhostBox(fluxSumE_(core::Component::X), [&](auto const&... args) mutable {
            if (std::isnan(fluxSumE_(core::Component::X)(args...)))
            {
                fluxSumE_(core::Component::X)(args...) = 0.0;
            }
            fluxSumE_(core::Component::X)(args...)
                += timeElectric(core::Component::X)(args...) * coef;
        });

        layout.evalOnGhostBox(fluxSumE_(core::Component::Y), [&](auto const&... args) mutable {
            if (std::isnan(fluxSumE_(core::Component::Y)(args...)))
            {
                fluxSumE_(core::Component::Y)(args...) = 0.0;
            }
            fluxSumE_(core::Component::Y)(args...)
                += timeElectric(core::Component::Y)(args...) * coef;
        });

        layout.evalOnGhostBox(fluxSumE_(core::Component::Z), [&](auto const&... args) mutable {
            if (std::isnan(fluxSumE_(core::Component::Z)(args...)))
            {
                fluxSumE_(core::Component::Z)(args...) = 0.0;
            }
            fluxSumE_(core::Component::Z)(args...)
                += timeElectric(core::Component::Z)(args...) * coef;
        });
    }
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger, ModelViews_t>::resetFluxSum(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level)
{
    auto& mhdModel = dynamic_cast<MHDModel&>(model);

    for (auto& patch : level)
    {
        auto const& layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto _             = mhdModel.resourcesManager->setOnPatch(*patch, fluxSum_, fluxSumE_);

        evalFluxesOnGhostBox(
            layout, [&](auto& left, auto const&... args) mutable { left(args...) = 0.0; },
            fluxSum_);

        layout.evalOnGhostBox(fluxSumE_(core::Component::X), [&](auto const&... args) mutable {
            fluxSumE_(core::Component::X)(args...) = 0.0;
        });

        layout.evalOnGhostBox(fluxSumE_(core::Component::Y), [&](auto const&... args) mutable {
            fluxSumE_(core::Component::Y)(args...) = 0.0;
        });

        layout.evalOnGhostBox(fluxSumE_(core::Component::Z), [&](auto const&... args) mutable {
            fluxSumE_(core::Component::Z)(args...) = 0.0;
        });
    }
}


template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger, ModelViews_t>::reflux(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, IMessenger& messenger,
    double const time)
{
    auto& fromCoarser                 = dynamic_cast<Messenger&>(messenger);
    auto& mhdModel                    = dynamic_cast<MHDModel&>(model);
    auto&& [timeFluxes, timeElectric] = evolve_.exposeFluxes();

    for (auto& patch : level)
    {
        core::FiniteVolumeEuler<GridLayout> fveuler;
        core::Faraday<GridLayout> faraday;

        auto& Bold_ = stateOld_.B;
        auto& B     = mhdModel.state.B;

        auto layout      = amr::layoutFromPatch<GridLayout>(*patch);
        auto _sp         = mhdModel.resourcesManager->setOnPatch(*patch, mhdModel.state, stateOld_,
                                                                 timeFluxes, Bold_, timeElectric, B);
        auto _sl_faraday = core::SetLayout(&layout, faraday);
        auto _sl_fveuler = core::SetLayout(&layout, fveuler);

        auto dt = time - oldTime_[level.getLevelNumber()];

        fveuler(stateOld_, mhdModel.state, timeFluxes, dt);
        faraday(Bold_, timeElectric, B, dt);
    };
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger,
         typename ModelViews_t>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger, ModelViews_t>::advanceLevel(
    hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& view,
    IMessenger& fromCoarserMessenger, double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::advanceLevel");

    auto& modelView   = dynamic_cast<ModelViews_t&>(view);
    auto& fromCoarser = dynamic_cast<Messenger&>(fromCoarserMessenger);
    auto level        = hierarchy.getPatchLevel(levelNumber);

    evolve_(modelView.model(), modelView.model().state, fluxes_, fromCoarser, *level, currentTime,
            newTime);
}

} // namespace PHARE::solver

#endif
