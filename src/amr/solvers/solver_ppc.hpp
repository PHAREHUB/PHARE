#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/algorithm.hpp"
#include "mpi/mpi_utils.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/physical_models/models.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/solvers/solver_hybrid_field_evolvers.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"

#include <SAMRAI/hier/Patch.h>

#include <memory>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <cmath>


namespace PHARE::solver
{
// -----------------------------------------------------------------------------


template<typename GridLayout, typename Boxing_t>
auto make_selection_boxes_for(auto const& hierarchy, auto& level)
{
    std::unordered_map<std::string, Boxing_t> levelBoxing;
    auto const neighbors = amr::makeSameLevelNeighbors(hierarchy, level.getLevelNumber());
    for (auto const& patch : level)
        if (auto [it, suc] = levelBoxing.try_emplace(
                core::to_string(patch->getGlobalId()),
                Boxing_t{amr::layoutFromPatch<GridLayout>(*patch),
                         amr::makeNonLevelGhostBoxFor<GridLayout>(*patch, neighbors),
                         amr::makeNonOwnedBoxesFor<GridLayout::dimension>(*patch, neighbors)});
            !suc)
            throw std::runtime_error("boxing map insertion failure");

    return levelBoxing;
}


// emits predictor1_/predictor2_/corrector_ in another translation unit, see solver_ppc_steps.hpp
template<typename Solver>
struct SolverPPCSteps;


template<typename HybridModel, typename AMR_Types>
class SolverPPC : public ISolver<AMR_Types>
{
    template<typename>
    friend struct SolverPPCSteps;

private:
    static constexpr auto dimension = HybridModel::dimension;

    using Electromag       = HybridModel::electromag_type;
    using Ions             = HybridModel::ions_type;
    using ParticleArray    = Ions::particle_array_type;
    using VecFieldT        = HybridModel::vecfield_type;
    using GridLayout       = HybridModel::gridlayout_type;
    using ResourcesManager = HybridModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using HybridMessenger  = amr::HybridMessenger<HybridModel>;

    using FE_t         = FieldEvolverDispatchers<HybridModel>;
    using Faraday_t    = FE_t::Faraday_t;
    using Ampere_t     = FE_t::Ampere_t;
    using Ohm_t        = OhmLevelTransformer<HybridModel>;
    using IonUpdater_t = core::IonUpdaterFor_t<ParticleArray, GridLayout>;
    using Boxing_t     = IonUpdater_t::Boxing_t;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

    VecFieldT Bold_{this->name() + "_Bold", core::HybridQuantity::Vector::B};
    VecFieldT fluxSumE_{this->name() + "_fluxSumE", core::HybridQuantity::Vector::E};
    std::unordered_map<std::size_t, double> oldTime_;

    core::OhmInfo ohm_info;
    IonUpdater_t ionUpdater_;


public:
    using patch_t     = AMR_Types::patch_t;
    using level_t     = AMR_Types::level_t;
    using hierarchy_t = AMR_Types::hierarchy_t;
    using level_ptr   = AMR_Types::level_ptr;


    explicit SolverPPC(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"PPC"}
        , ohm_info{core::OhmInfo::FROM(dict["ohm"])}
        , ionUpdater_{}

    {
    }

    ~SolverPPC() override = default;


    std::string modelName() const override { return HybridModel::model_name; }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    void registerResources(IPhysicalModel_t& model) override;


    void allocate(IPhysicalModel_t& model, patch_t& patch,
                  double const allocateTime) const override;

    void prepareStep(IPhysicalModel_t& model, level_t& level, double const currentTime) override;

    void accumulateFluxSum(IPhysicalModel_t& model, level_t& level, double const coef) override;

    void resetFluxSum(IPhysicalModel_t& model, level_t& level) override;

    void reflux(IPhysicalModel_t& model, level_t& level, IMessenger& messenger,
                double const time) override;

    double computeStableDt(IPhysicalModel_t& model, level_t& level,
                           CFLNumbers const& cflNumbers) override;

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, IPhysicalModel_t& model,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;


    void onRegrid() override
    {
        boxing.clear();
        ionUpdater_.reset();
    }


    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(Bold_, fluxSumE_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(Bold_, fluxSumE_);
    }


private:
    using Messenger = amr::HybridMessenger<HybridModel>;


    void predictor1_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void predictor2_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void corrector_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                    double const currentTime, double const newTime);


    void average_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                  double const currentTime, double const newTime);


    void moveIons_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                   double const currentTime, double const newTime, core::UpdaterMode mode);


    void make_boxes(hierarchy_t const& hierarchy, level_ptr const& level)
    {
        boxing[level->getLevelNumber()]
            = LevelBoxing{level, make_selection_boxes_for<GridLayout, Boxing_t>(hierarchy, *level)};
    }

    auto& setup_level(hierarchy_t const& hierarchy, int const levelNumber)
    {
        auto level = hierarchy.getPatchLevel(levelNumber);
        // level may have been removed and re-created without a regrid notification
        auto const it            = boxing.find(levelNumber);
        bool const needsCreating = it == boxing.end() or it->second.level.lock() != level;
        if (needsCreating)
            make_boxes(hierarchy, level);
        return *level;
    }

    void update_electrons(auto& level, auto& model)
    {
        ElectronsLevelTransformer{level, model}(model.state.electrons);
    }


    struct LevelBoxing
    {
        std::weak_ptr<level_t> level; // identity of the level the boxes were built for
        std::unordered_map<std::string /*patchid*/, Boxing_t> patches;
    };
    std::unordered_map<int /*level*/, LevelBoxing> boxing;


}; // end solverPPC



// -----------------------------------------------------------------------------



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::registerResources(IPhysicalModel_t& model)
{
    auto& hmodel = dynamic_cast<HybridModel&>(model);
    hmodel.resourcesManager->registerResources(electromagPred_);
    hmodel.resourcesManager->registerResources(electromagAvg_);

    hmodel.resourcesManager->registerResources(Bold_);
    hmodel.resourcesManager->registerResources(fluxSumE_);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::allocate(IPhysicalModel_t& model, patch_t& patch,
                                                 double const allocateTime) const
{
    auto& hmodel = dynamic_cast<HybridModel&>(model);
    hmodel.resourcesManager->allocate(electromagPred_, patch, allocateTime);
    hmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);

    hmodel.resourcesManager->allocate(Bold_, patch, allocateTime);
    hmodel.resourcesManager->allocate(fluxSumE_, patch, allocateTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    auto const& Eavg  = electromagAvg_.E;
    auto const& Bpred = electromagPred_.B;

    hybridInfo.ghostElectric.emplace_back(Eavg.name());
    hybridInfo.initMagnetic.emplace_back(Bpred.name());
    hybridInfo.ghostMagnetic.emplace_back(Bpred.name());
    hybridInfo.refluxElectric  = Eavg.name();
    hybridInfo.fluxSumElectric = fluxSumE_.name();
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::prepareStep(IPhysicalModel_t& model, level_t& level,
                                                    double const currentTime)
{
    oldTime_[level.getLevelNumber()] = currentTime;

    auto& hybridModel = dynamic_cast<HybridModel&>(model);

    for (auto& patch : level) // SAMRAI patch data times, serial
        hybridModel.resourcesManager->setTime(Bold_, *patch, currentTime);

    auto accessor
        = amr::make_model_level_accessor(level, hybridModel, hybridModel.state.electromag.B, Bold_);
    amr::level_exec(
        accessor, [](auto& view) -> auto& { return std::get<1>(view.args); },
        [](auto& view, auto const&... tile_idx) {
            auto& [B, Bold] = view.args;
            core::copy_data(Bold, B, tile_idx...);
        });
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::accumulateFluxSum(IPhysicalModel_t& model, level_t& level,
                                                          double const coef)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::accumulateFluxSum");

    auto& hybridModel = dynamic_cast<HybridModel&>(model);
    auto accessor = amr::make_model_level_accessor(level, hybridModel, fluxSumE_, electromagAvg_.E);
    amr::level_exec(
        accessor, [](auto& view) -> auto& { return std::get<0>(view.args); },
        [&](auto& view, auto const&... tile_idx) {
            auto& [fluxSumE, Eavg] = view.args;
            if constexpr (sizeof...(tile_idx) > 0)
            {
                auto dst       = core::tile_at(fluxSumE, tile_idx...);
                auto const src = core::tile_at(Eavg, tile_idx...);
                for (std::size_t c = 0; c < dst.size(); ++c)
                    core::operate_field<core::PlusEqualsProduct>(dst[c], src[c], coef);
            }
            else
                core::operate<core::PlusEqualsProduct>(fluxSumE, Eavg, coef);
        });
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::resetFluxSum(IPhysicalModel_t& model, level_t& level)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::resetFluxSum");

    auto& hybridModel = dynamic_cast<HybridModel&>(model);
    auto accessor     = amr::make_model_level_accessor(level, hybridModel, fluxSumE_);
    amr::level_exec(
        accessor, [](auto& view) -> auto& { return std::get<0>(view.args); },
        [](auto& view, auto const&... tile_idx) {
            core::maybe_tile_at(std::get<0>(view.args), tile_idx...).zero();
        });
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::reflux(IPhysicalModel_t& model, level_t& level,
                                               IMessenger& messenger, double const time)
{
    auto& hybridModel     = dynamic_cast<HybridModel&>(model);
    auto& hybridMessenger = dynamic_cast<HybridMessenger&>(messenger);
    auto& Eavg            = electromagAvg_.E;
    auto& B               = hybridModel.state.electromag.B;
    auto dt               = time - oldTime_[level.getLevelNumber()];

    Faraday_t{level, hybridModel}(Bold_, Eavg, B, dt);

    hybridMessenger.fillMagneticGhosts(B, level, time);
}


template<typename HybridModel, typename AMR_Types>
double SolverPPC<HybridModel, AMR_Types>::computeStableDt(IPhysicalModel_t& model, level_t& level,
                                                          CFLNumbers const& cflNumbers)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::computeStableDt");

    throw std::runtime_error("Adaptive time-stepping not ready in hybrid.");

    double const wave = cflNumbers.wave;

    auto& hybridModel = dynamic_cast<HybridModel&>(model);
    auto& n           = hybridModel.state.ions.massDensity();
    auto& B           = hybridModel.state.electromag.B;

    // whistler wave cfl criterion
    double dt = std::numeric_limits<double>::max();

    // on one patch, or one tile of it - tiles are traversed sequentially, so reducing into dt is
    //  safe; min is order independent and tile overlaps don't matter
    auto const whistlerDt = [&](auto const& layout, auto const& n, auto const& B) {
        auto const meshSize = layout.meshSize();

        auto const& Bx = B(core::Component::X);
        auto const& By = B(core::Component::Y);
        auto const& Bz = B(core::Component::Z);

        layout.evalOnBox(n, [&](auto&... args) mutable {
            core::MeshIndex<dimension> const index{args...};

            auto const rho = n(index);
            auto const bx  = GridLayout::template project<GridLayout::BxToMoments>(Bx, index);
            auto const by  = GridLayout::template project<GridLayout::ByToMoments>(By, index);
            auto const bz  = GridLayout::template project<GridLayout::BzToMoments>(Bz, index);

            auto const BdotB = bx * bx + by * by + bz * bz;

            double invDtWhistler = 0;
            for (std::size_t d = 0; d < dimension; ++d)
                invDtWhistler += std::numbers::pi
                                 * core::compute_whistler_(1.0 / meshSize[d], rho, BdotB)
                                 / meshSize[d];

            dt = std::min(dt, wave / invDtWhistler);
        });
    };

    auto& rm = *hybridModel.resourcesManager;
    for (auto& patch : rm.enumerate(level, n, B))
    {
        using field_type = std::decay_t<decltype(n)>;

        if constexpr (core::is_field_tile_set_v<field_type>)
            core::tile_exec_with_layout(whistlerDt, n, B);
        else
            whistlerDt(amr::layoutFromPatch<GridLayout>(*patch), n, B);
    }

    return dt;
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::advanceLevel(hierarchy_t const& hierarchy,
                                                     int const levelNumber,
                                                     IPhysicalModel_t& level_model,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::advanceLevel");

    auto& model       = dynamic_cast<HybridModel&>(level_model);
    auto& fromCoarser = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto& level       = setup_level(hierarchy, levelNumber);

    predictor1_(level, model, fromCoarser, currentTime, newTime);

    average_(level, model, fromCoarser, currentTime, newTime);

    moveIons_(level, model, fromCoarser, currentTime, newTime, core::UpdaterMode::domain_only);

    predictor2_(level, model, fromCoarser, currentTime, newTime);


    average_(level, model, fromCoarser, currentTime, newTime);

    moveIons_(level, model, fromCoarser, currentTime, newTime, core::UpdaterMode::all);

    corrector_(level, model, fromCoarser, currentTime, newTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, HybridModel& model,
                                                 Messenger& fromCoarser, double const currentTime,
                                                 double const newTime)
{
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::average_");

        TimeSetter setTime{level, model, newTime};
        auto accessor = amr::make_model_level_accessor(level, model, model.state.electromag,
                                                       electromagPred_, electromagAvg_);
        amr::level_exec(
            accessor, [](auto& view) -> auto& { return std::get<2>(view.args).B; },
            [](auto& view, auto const&... tile_idx) {
                auto& [electromag, pred, avg] = view.args;
                auto const average            = [&](auto& x, auto& x_pred, auto& x_avg) {
                    auto&& x_avg_t = core::maybe_tile_at(x_avg, tile_idx...);
                    PHARE::core::average(core::maybe_tile_at(x, tile_idx...),
                                                    core::maybe_tile_at(x_pred, tile_idx...), x_avg_t);
                };
                average(electromag.B, pred.B, avg.B);
                average(electromag.E, pred.E, avg.E);
            });

        setTime(electromagAvg_.B);
        setTime(electromagAvg_.E);
    }

    // the following will fill E on all edges of all ghost cells, including those
    // on domain border. For level ghosts, electric field will be obtained from
    // next coarser level E average

    PHARE_LOG_SCOPE(1, "SolverPPC::average::schedules");
    fromCoarser.fillElectricGhosts(electromagAvg_.E, level, (currentTime + newTime) / 2.);
}


void _debug_log_move_ions(auto& level, auto& model)
{
    auto& rm                              = *model.resourcesManager;
    std::size_t nbrLevelGhostOldParticles = 0;
    std::size_t nbrLevelGhostParticles    = 0;

    for (auto& patch : rm.enumerate(level, model.state.ions))
    {
        for (auto& pop : model.state.ions)
        {
            nbrLevelGhostOldParticles += pop.levelGhostParticlesOld().size();
            nbrLevelGhostParticles += pop.levelGhostParticles().size();

            if (nbrLevelGhostOldParticles < nbrLevelGhostParticles
                and nbrLevelGhostOldParticles > 0)
                throw std::runtime_error("Error - there are less old level ghost particles ("
                                         + std::to_string(nbrLevelGhostOldParticles)
                                         + ") than pushable ("
                                         + std::to_string(nbrLevelGhostParticles) + ")");
        }
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, HybridModel& model,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    using enum core::LayoutMode;
    using ParticleArray_t = Ions::particle_array_type;

    PHARE_DEBUG_DO(_debug_log_move_ions(level, model);)

    TimeSetter setTime{level, model, newTime};
    auto& rm                = *model.resourcesManager;
    auto const& levelBoxing = boxing[level.getLevelNumber()].patches;

    auto& ions = model.state.ions;
    try
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::moveIons");
        auto const dt = newTime - currentTime;
        if constexpr (ParticleArray_t::layout_mode == AoSMapped)
            for (auto& patch : rm.enumerate(level, ions, electromagAvg_))
                ionUpdater_.updatePopulations(ions, electromagAvg_,
                                              levelBoxing.at(core::to_string(patch->getGlobalId())),
                                              dt, mode);

        else
        {
            auto accessor = amr::make_model_level_accessor(level, model, ions, electromagAvg_);
            ionUpdater_.updatePopulations(accessor, levelBoxing, dt, mode);
        }
    }
    catch (core::DictionaryException const& ex)
    {
        PHARE_LOG_ERROR(ex());
    }
    if (mpi::any_errors())
        throw core::DictionaryException{}("ID", "Updater::updatePopulations");

    // this needs to be done before calling the messenger
    setTime(ions);
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::moveIons::fillFluxBorders");
        fromCoarser.fillFluxBorders(ions, level, newTime);
    }
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::moveIons::fillDensityBorders");
        fromCoarser.fillDensityBorders(ions, level, newTime);
    }
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::moveIons::fillIonPopMomentGhosts");
        fromCoarser.fillIonPopMomentGhosts(ions, level, newTime);
    }

    if (mode != core::UpdaterMode::domain_only)
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::moveIons::fillIonGhostParticles");

        // the filter below is exact only if the updater kept owned domain particles only
        PHARE_DEBUG_DO({
            for (auto& patch : rm.enumerate(level, ions))
            {
                auto const& boxes = levelBoxing.at(core::to_string(patch->getGlobalId()));
                for (auto& pop : ions)
                    core::per_particle(pop.domainParticles(), [&](auto const& particle) {
                        if (!boxes.isOwned(particle.iCell()))
                            throw std::runtime_error("domain particle outside owned cells");
                    });
            }
        })

        fromCoarser.fillIonGhostParticles(ions, level, newTime);

        // the exchange hands a leaving particle to every patch whose domain contains its
        // cell; on overlapping patches only the owner keeps it
        for (auto& patch : rm.enumerate(level, ions))
            ions.eraseDomainParticlesIn(
                levelBoxing.at(core::to_string(patch->getGlobalId())).nonOwnedBoxes);
    }

    {
        auto accessor = amr::make_model_level_accessor(level, model, ions);
        amr::level_exec(
            accessor, [](auto& view) -> auto& { return std::get<0>(view.args).chargeDensity(); },
            [](auto& view, auto const&... tile_idx) {
                std::get<0>(view.args).update(tile_idx...);
            });
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::moveIons::fillIonBorders");
        fromCoarser.fillIonBorders(ions, level, newTime);
    }
    // no need to update time, since it has been done before
}

} // namespace PHARE::solver

#endif
