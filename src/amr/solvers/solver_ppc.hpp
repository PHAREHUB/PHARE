#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_ppc_model_view.hpp"
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

    using ModelViews_t = HybridPPCModelView<HybridModel>;
    using Faraday_t    = ModelViews_t::Faraday_t;
    using Ampere_t     = ModelViews_t::Ampere_t;
    using Ohm_t        = ModelViews_t::Ohm_t;
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

    void reflux(IPhysicalModel_t& model, SAMRAI::hier::PatchLevel& level,
                double const time) override;

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& views,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;


    void onRegrid() override
    {
        boxing.clear();
        ionUpdater_.reset();
    }


    std::shared_ptr<ISolverModelView> make_view(level_t& level, IPhysicalModel_t& model) override
    {
        return std::make_shared<ModelViews_t>(level, dynamic_cast<HybridModel&>(model));
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


    void predictor1_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void predictor2_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void corrector_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                    double const currentTime, double const newTime);


    void average_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                  double const newTime);


    void moveIons_(level_t& level, ModelViews_t& views, Messenger& fromCoarser,
                   double const currentTime, double const newTime, core::UpdaterMode mode);


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
                                               SAMRAI::hier::PatchLevel& level, double const time)
{
    auto& hybridModel = dynamic_cast<HybridModel&>(model);
    auto& Eavg        = electromagAvg_.E;
    auto& B           = hybridModel.state.electromag.B;

    for (auto& patch : level)
    {
        core::Faraday<GridLayout> faraday;
        auto layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto _sp    = hybridModel.resourcesManager->setOnPatch(*patch, Bold_, Eavg, B);
        auto _sl    = core::SetLayout(&layout, faraday);
        auto dt     = time - oldTime_[level.getLevelNumber()];
        faraday(Bold_, Eavg, B, dt);
    };
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::advanceLevel(hierarchy_t const& hierarchy,
                                                     int const levelNumber, ISolverModelView& views,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::advanceLevel");

    auto& modelView   = dynamic_cast<ModelViews_t&>(views);
    auto& fromCoarser = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto& level       = setup_level(hierarchy, levelNumber);

    predictor1_(level, modelView, fromCoarser, currentTime, newTime);

    average_(level, modelView, fromCoarser, newTime);

    moveIons_(level, modelView, fromCoarser, currentTime, newTime, core::UpdaterMode::domain_only);

    predictor2_(level, modelView, fromCoarser, currentTime, newTime);


    average_(level, modelView, fromCoarser, newTime);

    moveIons_(level, modelView, fromCoarser, currentTime, newTime, core::UpdaterMode::all);

    corrector_(level, modelView, fromCoarser, currentTime, newTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, ModelViews_t& views,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_");

    TimeSetter setTime{views, newTime};

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.faraday");
        auto dt = newTime - currentTime;
        faraday_(views.layouts, views.electromag_B, views.electromag_E, views.electromagPred_B, dt);
        setTime([](auto& state) -> auto& { return state.electromagPred.B; });
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ampere");
        ampere_(views.layouts, views.electromagPred_B, views.J);
        setTime([](auto& state) -> auto& { return state.J; });
        fromCoarser.fillCurrentGhosts(views.model().state.J, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ohm");
        for (auto& state : views)
            state.electrons.update(state.layout);
        ohm_(views.layouts, views.N, views.Ve, views.Pe, views.electromagPred_B, views.J,
             views.electromagPred_E);
        setTime([](auto& state) -> auto& { return state.electromagPred.E; });
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, ModelViews_t& views,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_");

    TimeSetter setTime{views, newTime};

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.faraday");
        auto dt = newTime - currentTime;
        faraday_(views.layouts, views.electromag_B, views.electromagAvg_E, views.electromagPred_B,
                 dt);
        setTime([](auto& state) -> auto& { return state.electromagPred.B; });
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ampere");
        ampere_(views.layouts, views.electromagPred_B, views.J);
        setTime([](auto& state) -> auto& { return state.J; });
        fromCoarser.fillCurrentGhosts(views.model().state.J, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ohm");
        for (auto& state : views)
            state.electrons.update(state.layout);
        ohm_(views.layouts, views.N, views.Ve, views.Pe, views.electromagPred_B, views.J,
             views.electromagPred_E);
        setTime([](auto& state) -> auto& { return state.electromagPred.E; });
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, ModelViews_t& views,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::corrector_");

    auto levelNumber = level.getLevelNumber();
    TimeSetter setTime{views, newTime};

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.faraday");
        auto dt = newTime - currentTime;
        faraday_(views.layouts, views.electromag_B, views.electromagAvg_E, views.electromag_B, dt);
        setTime([](auto& state) -> auto& { return state.electromag.B; });
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ampere");
        ampere_(views.layouts, views.electromag_B, views.J);
        setTime([](auto& state) -> auto& { return state.J; });
        fromCoarser.fillCurrentGhosts(views.model().state.J, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm");
        for (auto& state : views)
            state.electrons.update(state.layout);
        ohm_(views.layouts, views.N, views.Ve, views.Pe, views.electromag_B, views.J,
             views.electromag_E);
        setTime([](auto& state) -> auto& { return state.electromag.E; });

        fromCoarser.fillElectricGhosts(views.model().state.electromag.E, level, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, ModelViews_t& views,
                                                 Messenger& fromCoarser, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::average_");

    TimeSetter setTime{views, newTime};

    for (auto& state : views)
    {
        PHARE::core::average(state.electromag.B, state.electromagPred.B, state.electromagAvg.B);
        PHARE::core::average(state.electromag.E, state.electromagPred.E, state.electromagAvg.E);
    }

    setTime([](auto& state) -> auto& { return state.electromagAvg.B; });
    setTime([](auto& state) -> auto& { return state.electromagAvg.E; });

    // the following will fill E on all edges of all ghost cells, including those
    // on domain border. For level ghosts, electric field will be obtained from
    // next coarser level E average
    fromCoarser.fillElectricGhosts(electromagAvg_.E, level, newTime);
}


template<typename... Args>
void _debug_log_move_ions(Args const&... args)
{
    auto const& [views] = std::forward_as_tuple(args...);

    std::size_t nbrLevelGhostOldParticles = 0;
    std::size_t nbrLevelGhostParticles    = 0;
    for (auto& state : views)
    {
        for (auto& pop : state.ions)
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
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, ModelViews_t& views,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::moveIons_");
    PHARE_DEBUG_DO(_debug_log_move_ions(views);)

    TimeSetter setTime{views, newTime};
    auto const& levelBoxing = boxing[level.getLevelNumber()];

    {
        auto dt = newTime - currentTime;
        for (auto& state : views)
            ionUpdater_.updatePopulations(
                state.ions, state.electromagAvg,
                levelBoxing.at(amr::to_string(state.patch->getGlobalId())), dt, mode);
    }

    // this needs to be done before calling the messenger
    setTime([](auto& state) -> auto& { return state.ions; });

    fromCoarser.fillFluxBorders(views.model().state.ions, level, newTime);
    fromCoarser.fillDensityBorders(views.model().state.ions, level, newTime);
    // fromCoarser.fillIonPopMomentGhosts(views.model().state.ions, level, newTime);
    fromCoarser.fillIonGhostParticles(views.model().state.ions, level, newTime);

    for (auto& state : views)
        ionUpdater_.updateIons(state.ions);
    // no need to update time, since it has been done before

    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    fromCoarser.fillIonMomentGhosts(views.model().state.ions, level, newTime);
}



} // namespace PHARE::solver




#endif
