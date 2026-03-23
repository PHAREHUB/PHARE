#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include "core/data/electrons/electrons.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"

#include <SAMRAI/hier/Patch.h>
#include "SAMRAI/hier/PatchLevel.h"

#include <tuple>
#include <unordered_map>


namespace PHARE::solver
{
// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
class SolverPPC : public ISolver<AMR_Types>
{
private:
    static constexpr auto dimension    = HybridModel::dimension;
    static constexpr auto interp_order = HybridModel::gridlayout_type::interp_order;

    using Electromag       = HybridModel::electromag_type;
    using Ions             = HybridModel::ions_type;
    using ParticleArray    = Ions::particle_array_type;
    using VecFieldT        = HybridModel::vecfield_type;
    using GridLayout       = HybridModel::gridlayout_type;
    using ResourcesManager = HybridModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using HybridMessenger  = amr::HybridMessenger<HybridModel>;

    using Faraday_t    = FaradayTransformer<GridLayout, AMR_Types>;
    using Ampere_t     = AmpereTransformer<GridLayout, AMR_Types>;
    using Ohm_t        = OhmTransformer<GridLayout, AMR_Types>;
    using IonUpdater_t = PHARE::core::IonUpdater<Ions, Electromag, GridLayout>;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

    VecFieldT Bold_{this->name() + "_Bold", core::HybridQuantity::Vector::B};
    VecFieldT fluxSumE_{this->name() + "_fluxSumE", core::HybridQuantity::Vector::E};
    std::unordered_map<std::size_t, double> oldTime_;

    Faraday_t faraday_;
    Ampere_t ampere_;
    Ohm_t ohm_;

    IonUpdater_t ionUpdater_;


public:
    using patch_t     = AMR_Types::patch_t;
    using level_t     = AMR_Types::level_t;
    using hierarchy_t = AMR_Types::hierarchy_t;



    explicit SolverPPC(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"PPC"}
        , ohm_{dict["ohm"]}
        , ionUpdater_{dict["ion_updater"]}

    {
    }

    ~SolverPPC() override = default;


    std::string modelName() const override { return HybridModel::model_name; }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    void registerResources(IPhysicalModel_t& model) override;


    void allocate(IPhysicalModel_t& model, SAMRAI::hier::Patch& patch,
                  double const allocateTime) const override;

    void prepareStep(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                     double const currentTime) override;

    void accumulateFluxSum(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                           double const coef) override;

    void resetFluxSum(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level) override;

    void reflux(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level, IMessenger& messenger,
                double const time) override;

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, IPhysicalModel_t& views,
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


    void average_(level_t& level, HybridModel& model, Messenger& fromCoarser, double const newTime);


    void moveIons_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                   double const currentTime, double const newTime, core::UpdaterMode mode);


    struct TimeSetter
    {
        void operator()(auto&... quantities)
        {
            auto& rm = *model.resourcesManager;
            for (auto& patch : rm.enumerate(level, quantities...))
                (model.resourcesManager->setTime(quantities, *patch, newTime), ...);
        }

        level_t& level;
        HybridModel& model;
        double newTime;
    };


    void make_boxes(hierarchy_t const& hierarchy, level_t& level)
    {
        int const lvlNbr = level.getLevelNumber();
        if (boxing.count(lvlNbr))
            return;

        auto& levelBoxing = boxing[lvlNbr]; // creates if missing

        for (auto const& patch : level)
            if (auto [it, suc] = levelBoxing.try_emplace(
                    amr::to_string(patch->getGlobalId()),
                    Boxing_t{amr::layoutFromPatch<GridLayout>(*patch),
                             amr::makeNonLevelGhostBoxFor<GridLayout>(*patch, hierarchy)});
                !suc)
                throw std::runtime_error("boxing map insertion failure");
    }

    auto& setup_level(hierarchy_t const& hierarchy, int const levelNumber)
    {
        auto level = hierarchy.getPatchLevel(levelNumber);
        if (boxing.count(levelNumber) == 0)
            make_boxes(hierarchy, *level);
        return *level;
    }

    void update_electrons(auto& level, auto& model)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, model.state.electrons))
            model.state.electrons.update(amr::layoutFromPatch<GridLayout>(*patch));
    };


    using Boxing_t = core::UpdaterSelectionBoxing<IonUpdater_t, GridLayout>;
    std::unordered_map<int /*level*/, std::unordered_map<std::string /*patchid*/, Boxing_t>> boxing;


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
void SolverPPC<HybridModel, AMR_Types>::allocate(IPhysicalModel_t& model,
                                                 SAMRAI::hier::Patch& patch,
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
void SolverPPC<HybridModel, AMR_Types>::prepareStep(IPhysicalModel_t& model,
                                                    SAMRAI::hier::PatchLevel& level,
                                                    double const currentTime)
{
    oldTime_[level.getLevelNumber()] = currentTime;

    auto& hybridModel = dynamic_cast<HybridModel&>(model);
    auto& B           = hybridModel.state.electromag.B;

    for (auto& patch : level)
    {
        auto dataOnPatch = hybridModel.resourcesManager->setOnPatch(*patch, B, Bold_);

        hybridModel.resourcesManager->setTime(Bold_, *patch, currentTime);

        Bold_.copyData(B);
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::accumulateFluxSum(IPhysicalModel_t& model,
                                                          SAMRAI::hier::PatchLevel& level,
                                                          double const coef)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::accumulateFluxSum");

    auto& hybridModel = dynamic_cast<HybridModel&>(model);

    for (auto& patch : level)
    {
        auto& Eavg         = electromagAvg_.E;
        auto const& layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto _             = hybridModel.resourcesManager->setOnPatch(*patch, fluxSumE_, Eavg);

        layout.evalOnGhostBox(fluxSumE_(core::Component::X), [&](auto const&... args) mutable {
            fluxSumE_(core::Component::X)(args...) += Eavg(core::Component::X)(args...) * coef;
        });

        layout.evalOnGhostBox(fluxSumE_(core::Component::Y), [&](auto const&... args) mutable {
            fluxSumE_(core::Component::Y)(args...) += Eavg(core::Component::Y)(args...) * coef;
        });

        layout.evalOnGhostBox(fluxSumE_(core::Component::Z), [&](auto const&... args) mutable {
            fluxSumE_(core::Component::Z)(args...) += Eavg(core::Component::Z)(args...) * coef;
        });
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::resetFluxSum(IPhysicalModel_t& model,
                                                     SAMRAI::hier::PatchLevel& level)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::resetFluxSum");

    auto& hybridModel = dynamic_cast<HybridModel&>(model);

    for (auto& patch : level)
    {
        auto const& layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto _             = hybridModel.resourcesManager->setOnPatch(*patch, fluxSumE_);

        fluxSumE_.zero();
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::reflux(IPhysicalModel_t& model,
                                               SAMRAI::hier::PatchLevel& level,
                                               IMessenger& messenger, double const time)
{
    auto& hybridModel     = dynamic_cast<HybridModel&>(model);
    auto& hybridMessenger = dynamic_cast<HybridMessenger&>(messenger);
    auto& Eavg            = electromagAvg_.E;
    auto& B               = hybridModel.state.electromag.B;
    auto dt               = time - oldTime_[level.getLevelNumber()];
    faraday_(level, hybridModel, Bold_, Eavg, B, dt);

    hybridMessenger.fillMagneticGhosts(B, level, time);
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

    average_(level, model, fromCoarser, newTime);

    moveIons_(level, model, fromCoarser, currentTime, newTime, core::UpdaterMode::domain_only);

    predictor2_(level, model, fromCoarser, currentTime, newTime);


    average_(level, model, fromCoarser, newTime);

    moveIons_(level, model, fromCoarser, currentTime, newTime, core::UpdaterMode::all);

    corrector_(level, model, fromCoarser, currentTime, newTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_");

    TimeSetter setTime{level, model, newTime};

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.faraday");
        auto dt = newTime - currentTime;
        faraday_(level, model, model.state.electromag.B, model.state.electromag.E,
                 electromagPred_.B, dt);
        setTime(electromagPred_.B);
        fromCoarser.fillMagneticGhosts(electromagPred_.B, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ampere");
        ampere_(level, model, electromagPred_.B, model.state.J);
        setTime(model.state.J);
        fromCoarser.fillCurrentGhosts(model.state.J, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ohm");
        update_electrons(level, model);
        ohm_(level, model, electromagPred_.B, electromagPred_.E);
        setTime(electromagPred_.E);
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_");

    TimeSetter setTime{level, model, newTime};

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.faraday");
        auto dt = newTime - currentTime;
        faraday_(level, model, model.state.electromag.B, electromagAvg_.E, electromagPred_.B, dt);
        setTime(electromagPred_.B);
        fromCoarser.fillMagneticGhosts(electromagPred_.B, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ampere");
        ampere_(level, model, electromagPred_.B, model.state.J);
        setTime(model.state.J);
        fromCoarser.fillCurrentGhosts(model.state.J, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ohm");
        update_electrons(level, model);
        ohm_(level, model, electromagPred_.B, electromagPred_.E);
        setTime(electromagPred_.E);
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, HybridModel& model,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::corrector_");

    auto levelNumber = level.getLevelNumber();
    TimeSetter setTime{level, model, newTime};

    auto& electromag = model.state.electromag;
    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.faraday");
        auto dt = newTime - currentTime;
        faraday_(level, model, electromag.B, electromagAvg_.E, electromag.B, dt);
        setTime(model.state.electromag.B);
        fromCoarser.fillMagneticGhosts(model.state.electromag.B, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ampere");
        ampere_(level, model, electromag.B, model.state.J);
        setTime(model.state.J);
        fromCoarser.fillCurrentGhosts(model.state.J, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm");

        update_electrons(level, model);
        ohm_(level, model, electromag.B, electromag.E);
        setTime(model.state.electromag.E);

        fromCoarser.fillElectricGhosts(model.state.electromag.E, level, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, HybridModel& model,
                                                 Messenger& fromCoarser, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::average_");

    TimeSetter setTime{level, model, newTime};
    auto& rm         = *model.resourcesManager;
    auto& electromag = model.state.electromag;
    for (auto& _ : rm.enumerate(level, electromag, electromagPred_, electromagAvg_))
    {
        PHARE::core::average(electromag.B, electromagPred_.B, electromagAvg_.B);
        PHARE::core::average(electromag.E, electromagPred_.E, electromagAvg_.E);
    }

    setTime(electromagAvg_.B);
    setTime(electromagAvg_.E);

    // the following will fill E on all edges of all ghost cells, including those
    // on domain border. For level ghosts, electric field will be obtained from
    // next coarser level E average
    fromCoarser.fillElectricGhosts(electromagAvg_.E, level, newTime);
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
    PHARE_LOG_SCOPE(1, "SolverPPC::moveIons_");
    PHARE_DEBUG_DO(_debug_log_move_ions(level, model);)

    TimeSetter setTime{level, model, newTime};
    auto& rm                = *model.resourcesManager;
    auto const& levelBoxing = boxing[level.getLevelNumber()];

    auto& ions = model.state.ions;
    try
    {
        auto dt = newTime - currentTime;
        for (auto& patch : rm.enumerate(level, ions, electromagAvg_))
            ionUpdater_.updatePopulations(ions, electromagAvg_,
                                          levelBoxing.at(amr::to_string(patch->getGlobalId())), dt,
                                          mode);
    }
    catch (core::DictionaryException const& ex)
    {
        PHARE_LOG_ERROR(ex());
    }
    if (core::mpi::any(core::Errors::instance().any()))
        throw core::DictionaryException{}("ID", "Updater::updatePopulations");

    // this needs to be done before calling the messenger
    setTime(ions);

    fromCoarser.fillFluxBorders(model.state.ions, level, newTime);
    fromCoarser.fillDensityBorders(model.state.ions, level, newTime);
    fromCoarser.fillIonPopMomentGhosts(model.state.ions, level, newTime);
    fromCoarser.fillIonGhostParticles(model.state.ions, level, newTime);

    for (auto& patch : rm.enumerate(level, ions))
        ionUpdater_.updateIons(ions);

    fromCoarser.fillIonBorders(model.state.ions, level, newTime);

    // no need to update time, since it has been done before
    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    // fromCoarser.fillIonMomentGhosts(model.state.ions, level, newTime);
}



} // namespace PHARE::solver




#endif
