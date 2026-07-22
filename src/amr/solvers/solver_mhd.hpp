#ifndef PHARE_SOLVER_MHD_HPP
#define PHARE_SOLVER_MHD_HPP


#include "core/def.hpp"
#include "initializer/data_provider.hpp"

#include "core/errors.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/mhd_messenger.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver_mhd_field_evolvers.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"


#include <array>
#include <cmath>
#include <tuple>
#include <vector>
#include <stdexcept>
#include <functional>
#include <type_traits>

namespace PHARE::solver
{
template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy,
         typename Messenger = amr::MHDMessenger<MHDModel>>
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
    EulerUsingComputedFlux<MHDModel> reflux_euler_;

    std::unordered_map<std::size_t, double> oldTime_;

    // adaptive-timestep coefficients (read from the algo dict, mirror ComputeFluxes/CT keys)
    double const gamma_; // adiabatic index (advective fast speed)
    double const eta_;   // resistivity (parabolic / Fourier bucket)
    bool const hall_;    // Hall active -> add whistler speed to the advective bucket

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
        , gamma_{dict["to_primitive"]["heat_capacity_ratio"].template to<double>()}
        , eta_{dict["constrained_transport"]["resistivity"].template to<double>()}
        , hall_{cppdict::get_value(dict, "fv_method/hall", false)}
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

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, IPhysicalModel_t& model,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;

    double computeStableDt(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                           StabilityNumbers const& stability) override;

    void onRegrid() override {}


    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(fluxes_, fluxSum_, fluxSumE_, stateOld_, evolve_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(fluxes_, fluxSum_, fluxSumE_, stateOld_, evolve_);
    }

private:
    void mhdNaNCheck_(MHDModel& state, level_t const& level, double time);
};

// -----------------------------------------------------------------------------

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::registerResources(
    IPhysicalModel_t& model)
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

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::allocate(
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

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const

{
    auto& mhdInfo = dynamic_cast<amr::MHDMessengerInfo&>(*info);

    mhdInfo.ghostMagneticFluxesX.emplace_back(fluxes_.B_fx.name());

    if constexpr (dimension >= 2)
    {
        mhdInfo.ghostMagneticFluxesY.emplace_back(fluxes_.B_fy.name());

        if constexpr (dimension == 3)
        {
            mhdInfo.ghostMagneticFluxesZ.emplace_back(fluxes_.B_fz.name());
        }
    }

    evolve_.fillMessengerInfo(mhdInfo);

    auto&& [timeFluxes, timeElectric] = evolve_.exposeFluxes();

    mhdInfo.reflux          = core::AllFluxesNames{timeFluxes};
    mhdInfo.refluxElectric  = timeElectric.name();
    mhdInfo.fluxSum         = core::AllFluxesNames{fluxSum_};
    mhdInfo.fluxSumElectric = fluxSumE_.name();

    // for the faraday in reflux
    mhdInfo.ghostElectric.emplace_back(timeElectric.name());
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::prepareStep(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, double const currentTime)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::prepareStep");

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


template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::accumulateFluxSum(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, double const coef)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::accumulateFluxSum");

    auto& mhdModel = dynamic_cast<MHDModel&>(model);
    auto& rm       = *mhdModel.resourcesManager;

    // MacOS clang still unhappy with structured bindings captures in lambdas
    auto&& tf          = evolve_.exposeFluxes();
    auto& timeFluxes   = std::get<0>(tf);
    auto& timeElectric = std::get<1>(tf);

    for (auto& _ : rm.enumerate(level, fluxSum_, fluxSumE_, timeFluxes, timeElectric))
    {
        core::operate<core::PlusEqualsProduct>(fluxSum_, timeFluxes, coef);
        core::operate<core::PlusEqualsProduct>(fluxSumE_, timeElectric, coef);
    }
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::resetFluxSum(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level)
{
    auto& mhdModel = dynamic_cast<MHDModel&>(model);
    auto& rm       = *mhdModel.resourcesManager;

    for (auto& _ : rm.enumerate(level, fluxSum_, fluxSumE_))
    {
        fluxSum_.zero();
        fluxSumE_.zero();
    }
}


template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::reflux(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, IMessenger& messenger,
    double const time)
{
    auto& bc                          = dynamic_cast<Messenger&>(messenger);
    auto& mhdModel                    = dynamic_cast<MHDModel&>(model);
    auto&& [timeFluxes, timeElectric] = evolve_.exposeFluxes();

    reflux_euler_(mhdModel, stateOld_, mhdModel.state, timeElectric, timeFluxes, bc, level, time,
                  time - oldTime_[level.getLevelNumber()]);
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::advanceLevel(
    hierarchy_t const& hierarchy, int const levelNumber, IPhysicalModel_t& model,
    IMessenger& fromCoarserMessenger, double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::advanceLevel");

    auto& mhdModel    = dynamic_cast<MHDModel&>(model);
    auto& fromCoarser = dynamic_cast<Messenger&>(fromCoarserMessenger);
    auto level        = hierarchy.getPatchLevel(levelNumber);

    try
    {
        evolve_(mhdModel, mhdModel.state, fluxes_, fromCoarser, *level, currentTime, newTime);

        PHARE_DEBUG_DO({ mhdNaNCheck_(mhdModel, *level, currentTime); })
    }
    catch (core::DictionaryException& ex)
    {
        PHARE_LOG_ERROR(ex());
    }

    if (core::mpi::any_errors())
        throw core::DictionaryException{}("ID", "SolverMHD::advanceLevel");
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
double SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::computeStableDt(
    IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, StabilityNumbers const& stability)
{
    PHARE_LOG_SCOPE(1, "SolverMHD::computeStableDt");

    auto const [cfl, fourier] = stability;

    auto& mhdModel = dynamic_cast<MHDModel&>(model);
    auto& rho      = mhdModel.state.rho;
    auto& rhoV     = mhdModel.state.rhoV;
    auto& B        = mhdModel.state.B;
    auto& Etot     = mhdModel.state.Etot;

    // Two stability buckets, combined by min. Both coefficients are normalized so that the value
    // 1 sits exactly on the (forward-Euler / SSP-RK) stability limit, independent of dimension, so
    // cfl, fourier are meant to be chosen in (0, 1]:
    //   - advective: dt = cfl / sum_d (|v_d| + c_fast_d [+ c_whistler_d if Hall]) / dx_d
    //   - resistive: dt = fourier / (2 * eta * sum_d 1/dx_d^2)   (eta uniform)
    // The level's patches are distributed across ranks, so the local min below is reduced across
    // ranks before returning. The inter-level projection is applied by the caller.
    double dt = std::numeric_limits<double>::max();

    for (auto& patch : level)
    {
        auto const& layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto _             = mhdModel.resourcesManager->setOnPatch(*patch, rho, rhoV, B, Etot);

        auto const meshSize = layout.meshSize();

        // resistive (Fourier) bucket: eta uniform -> one value per patch, no cell loop needed
        if (eta_ > 0)
        {
            double invdx2 = 0;
            for (std::size_t d = 0; d < dimension; ++d)
                invdx2 += 1.0 / (meshSize[d] * meshSize[d]);
            dt = std::min(dt, fourier / (2.0 * eta_ * invdx2));
        }

        auto const& rhoVx = rhoV(core::Component::X);
        auto const& rhoVy = rhoV(core::Component::Y);
        auto const& rhoVz = rhoV(core::Component::Z);
        auto const& Bx    = B(core::Component::X);
        auto const& By    = B(core::Component::Y);
        auto const& Bz    = B(core::Component::Z);

        // advective (+ Hall whistler) bucket: per cell, sum-of-speeds form
        layout.evalOnBox(rho, [&](auto&... args) mutable {
            core::MeshIndex<dimension> const index{args...};

            auto const r  = rho(index);
            auto const vx = rhoVx(index) / r;
            auto const vy = rhoVy(index) / r;
            auto const vz = rhoVz(index) / r;
            // cell-center the face-centered (Yee) fields, same idiom as ToPrimitiveConverter
            auto const bx = GridLayout::template project<GridLayout::faceXToCellCenter>(Bx, index);
            auto const by = GridLayout::template project<GridLayout::faceYToCellCenter>(By, index);
            auto const bz = GridLayout::template project<GridLayout::faceZToCellCenter>(Bz, index);

            auto const P     = core::eosEtotToP(gamma_, r, vx, vy, vz, bx, by, bz, Etot(index));
            auto const BdotB = bx * bx + by * by + bz * bz;

            std::array<double, 3> const v{vx, vy, vz};
            std::array<double, 3> const b{bx, by, bz};

            // sum_d (|v_d| + c_fast_d + c_whistler_d) / dx_d over simulated directions
            double invDtAdv = 0;
            for (std::size_t d = 0; d < dimension; ++d)
            {
                auto const cfast = core::compute_fast_magnetosonic_(gamma_, r, b[d], BdotB, P);
                // Hall whistler
                auto const cw = hall_ ? core::compute_whistler_(1.0 / meshSize[d], r, BdotB) : 0.0;
                invDtAdv += (std::abs(v[d]) + cfast + cw) / meshSize[d];
            }
            dt = std::min(dt, cfl / invDtAdv);
        });
    }

    return dt; // LOCAL (per-rank) min; caller reduces once across ranks for the whole cascade
}

template<typename MHDModel, typename AMR_Types, typename TimeIntegratorStrategy, typename Messenger>
void SolverMHD<MHDModel, AMR_Types, TimeIntegratorStrategy, Messenger>::mhdNaNCheck_(
    MHDModel& model, level_t const& level, double time)
{
    auto& rm  = model.resourcesManager;
    auto& rho = model.state.rho;

    auto check_nans = [&](auto const& field, auto const& origin,
                          core::MeshIndex<MHDModel::dimension> const& index) {
        if (std::isnan(field(index)))
        {
            std::stringstream ss;
            ss << "NaN detected in MHD field at index " << index << " on patch of origin " << origin
               << " on level " << level.getLevelNumber() << " at time " << time;
            core::DictionaryException ex{"cause", ss.str()};
            throw ex;
        }
    };

    for (auto const& patch : rm->enumerate(level, rho))
    {
        auto layout = amr::layoutFromPatch<GridLayout>(*patch);
        layout.evalOnGhostBox(
            rho, [&](auto const&... args) { check_nans(rho, layout.origin(), {args...}); });
    }
}

} // namespace PHARE::solver

#endif
