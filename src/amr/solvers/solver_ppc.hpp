#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/Patch.h>


#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_ppc_model_view.hpp"

#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/ohm/ohm.hpp"

#include "core/data/vecfield/vecfield.hpp"
#include "core/data/grid/gridlayout_utils.hpp"


#include <iomanip>
#include <sstream>

namespace PHARE::solver
{
// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
class SolverPPC : public ISolver<AMR_Types>
{
private:
    static constexpr auto dimension    = HybridModel::dimension;
    static constexpr auto interp_order = HybridModel::gridlayout_type::interp_order;

    using Electromag       = typename HybridModel::electromag_type;
    using Ions             = typename HybridModel::ions_type;
    using ParticleArray    = typename Ions::particle_array_type;
    using VecFieldT        = typename HybridModel::vecfield_type;
    using GridLayout       = typename HybridModel::gridlayout_type;
    using ResourcesManager = typename HybridModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using HybridMessenger  = amr::HybridMessenger<HybridModel>;

    using HybridModelView_t = HybridPPCModelView<HybridModel>;
    using Faraday_t         = typename HybridModelView_t::Faraday_t;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};


    Faraday_t faraday_1;
    PHARE::core::Faraday<GridLayout> faraday_;
    PHARE::core::Ampere<GridLayout> ampere_;
    PHARE::core::Ohm<GridLayout> ohm_;
    PHARE::core::IonUpdater<Ions, Electromag, GridLayout> ionUpdater_;



public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;



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



    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& view,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;

    std::shared_ptr<ISolverModelView> make_view(level_t& level, IPhysicalModel_t& model) override
    {
        return std::make_shared<HybridModelView_t>(level, dynamic_cast<HybridModel&>(model));
    }

private:
    using Messenger = amr::HybridMessenger<HybridModel>;


    void predictor1_(level_t& level, HybridModelView_t& view, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void predictor2_(level_t& level, HybridModelView_t& view, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void corrector_(level_t& level, HybridModelView_t& view, Messenger& fromCoarser,
                    double const currentTime, double const newTime);


    void average_(level_t& level, HybridModelView_t& view, Messenger& fromCoarser,
                  double const newTime);


    void moveIons_(level_t& level, HybridModelView_t& view, Messenger& fromCoarser,
                   double const currentTime, double const newTime, core::UpdaterMode mode);


    void saveState_(level_t& level, HybridModelView_t& view);

    void restoreState_(level_t& level, HybridModelView_t& view);


    // extend lifespan
    std::unordered_map<std::string, ParticleArray> tmpDomain;
    std::unordered_map<std::string, ParticleArray> patchGhost;


}; // end solverPPC



// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::registerResources(IPhysicalModel_t& model)
{
    auto& hmodel = dynamic_cast<HybridModel&>(model);
    hmodel.resourcesManager->registerResources(electromagPred_);
    hmodel.resourcesManager->registerResources(electromagAvg_);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::allocate(IPhysicalModel_t& model,
                                                 SAMRAI::hier::Patch& patch,
                                                 double const allocateTime) const
{
    auto& hmodel = dynamic_cast<HybridModel&>(model);
    hmodel.resourcesManager->allocate(electromagPred_, patch, allocateTime);
    hmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    auto const& Eavg  = electromagAvg_.E;
    auto const& Bpred = electromagPred_.B;

    hybridInfo.ghostElectric.emplace_back(core::VecFieldNames{Eavg});
    hybridInfo.initMagnetic.emplace_back(core::VecFieldNames{Bpred});
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::saveState_(level_t& level, HybridModelView_t& view)
{
    for (auto& state : view)
    {
        std::stringstream ss;
        ss << state.patch->getGlobalId();
        for (auto& pop : state.ions)
        {
            auto key = ss.str() + "_" + pop.name();
            if (!tmpDomain.count(key))
                tmpDomain.emplace(key, pop.domainParticles());
            else
                tmpDomain.at(key) = pop.domainParticles();
            if (!patchGhost.count(key))
                patchGhost.emplace(key, pop.patchGhostParticles());
            else
                patchGhost.at(key) = pop.patchGhostParticles();
        }
    }
}

template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::restoreState_(level_t& level, HybridModelView_t& view)
{
    for (auto& state : view)
    {
        std::stringstream ss;
        ss << state.patch->getGlobalId();

        for (auto& pop : state.ions)
        {
            pop.domainParticles()     = std::move(tmpDomain.at(ss.str() + "_" + pop.name()));
            pop.patchGhostParticles() = std::move(patchGhost.at(ss.str() + "_" + pop.name()));
        }
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::advanceLevel(hierarchy_t const& hierarchy,
                                                     int const levelNumber, ISolverModelView& view,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::advanceLevel");

    auto& modelView   = dynamic_cast<HybridModelView_t&>(view);
    auto& fromCoarser = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto level        = hierarchy.getPatchLevel(levelNumber);

    predictor1_(*level, modelView, fromCoarser, currentTime, newTime);

    average_(*level, modelView, fromCoarser, newTime);

    saveState_(*level, modelView);

    moveIons_(*level, modelView, fromCoarser, currentTime, newTime, core::UpdaterMode::domain_only);

    predictor2_(*level, modelView, fromCoarser, currentTime, newTime);


    average_(*level, modelView, fromCoarser, newTime);

    restoreState_(*level, modelView);

    moveIons_(*level, modelView, fromCoarser, currentTime, newTime, core::UpdaterMode::all);

    corrector_(*level, modelView, fromCoarser, currentTime, newTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, HybridModelView_t& view,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::predictor1_");

    auto dt                = newTime - currentTime;
    auto& resourcesManager = view.model().resourcesManager;

    {
        PHARE_LOG_SCOPE("SolverPPC::predictor1_.faraday");

        faraday_1(
            view,
            [](auto& state) {
                return std::forward_as_tuple( //
                    state.layout, state.electromag.B, state.electromag.E, state.electromagPred.B);
            },
            dt);

        for (auto& state : view)
            resourcesManager->setTime(state.electromagPred.B, *state.patch, newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::predictor1_.ampere");
        ampere_.op(view, [](auto& state) { return state.ppc_predictor1_ampere(); });
        for (auto& state : view)
            resourcesManager->setTime(state.J, *state.patch, newTime);
        fromCoarser.fillCurrentGhosts(view.model().state.J, level.getLevelNumber(), newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::predictor1_.ohm");
        for (auto& state : view)
            state.electrons.update(state.layout);
        ohm_(view, [](auto& state) { return state.ppc_predictor1_ohm(); });
        for (auto& state : view)
            resourcesManager->setTime(state.electromagPred.E, *state.patch, newTime);
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, HybridModelView_t& view,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::predictor2_");

    auto dt                = newTime - currentTime;
    auto& resourcesManager = view.model().resourcesManager;

    {
        PHARE_LOG_SCOPE("SolverPPC::predictor2_.faraday");
        faraday_(
            view, [](auto& state) { return state.ppc_predictor2_faraday(); }, dt);
        for (auto& state : view)
            resourcesManager->setTime(state.electromagPred.B, *state.patch, newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::predictor2_.ampere");
        ampere_.op(view, [](auto& state) { return state.ppc_predictor2_ampere(); });
        for (auto& state : view)
            resourcesManager->setTime(state.J, *state.patch, newTime);
        fromCoarser.fillCurrentGhosts(view.model().state.J, level.getLevelNumber(), newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::predictor2_.ohm");
        for (auto& state : view)
            state.electrons.update(state.layout);
        ohm_(view, [](auto& state) { return state.ppc_predictor2_ohm(); });
        for (auto& state : view)
            resourcesManager->setTime(state.electromagPred.E, *state.patch, newTime);
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, HybridModelView_t& view,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::corrector_");

    auto dt                = newTime - currentTime;
    auto levelNumber       = level.getLevelNumber();
    auto& resourcesManager = view.model().resourcesManager;

    {
        PHARE_LOG_SCOPE("SolverPPC::corrector_.faraday");
        faraday_(
            view, [](auto& state) { return state.ppc_corrector_faraday(); }, dt);
        for (auto& state : view)
            resourcesManager->setTime(state.electromag.B, *state.patch, newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::corrector_.ampere");
        ampere_.op(view, [](auto& state) { return state.ppc_corrector_ampere(); });
        for (auto& state : view)
            resourcesManager->setTime(state.J, *state.patch, newTime);
        fromCoarser.fillCurrentGhosts(view.model().state.J, levelNumber, newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::corrector_.ohm");
        for (auto& state : view)
            state.electrons.update(state.layout);
        ohm_(view, [](auto& state) { return state.ppc_corrector_ohm(); });
        for (auto& state : view)
            resourcesManager->setTime(state.electromag.E, *state.patch, newTime);

        fromCoarser.fillElectricGhosts(view.model().state.electromag.E, levelNumber, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, HybridModelView_t& view,
                                                 Messenger& fromCoarser, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::average_");

    auto _average = [](auto&... args) { PHARE::core::average(args...); };
    for (auto& state : view)
    {
        std::apply(_average, state.ppc_average_B());
        std::apply(_average, state.ppc_average_E());
    }

    // the following will fill E on all edges of all ghost cells, including those
    // on domain border. For level ghosts, electric field will be obtained from
    // next coarser level E average
    fromCoarser.fillElectricGhosts(electromagAvg_.E, level.getLevelNumber(), newTime);
}


template<typename... Args>
void _debug_log_move_ions(Args const&... args)
{
    auto const& [view] = std::forward_as_tuple(args...);

    std::size_t nbrDomainParticles        = 0;
    std::size_t nbrPatchGhostParticles    = 0;
    std::size_t nbrLevelGhostNewParticles = 0;
    std::size_t nbrLevelGhostOldParticles = 0;
    std::size_t nbrLevelGhostParticles    = 0; //
    for (auto& state : view)
    {
        for (auto& pop : state.ions)
        {
            nbrDomainParticles += pop.domainParticles().size();
            nbrPatchGhostParticles += pop.patchGhostParticles().size();
            nbrLevelGhostNewParticles += pop.levelGhostParticlesNew().size();
            nbrLevelGhostOldParticles += pop.levelGhostParticlesOld().size();
            nbrLevelGhostParticles += pop.levelGhostParticles().size();
            nbrPatchGhostParticles += pop.patchGhostParticles().size();

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
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, HybridModelView_t& view,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    PHARE_LOG_SCOPE("SolverPPC::moveIons_");

    PHARE_DEBUG_DO(_debug_log_move_ions(view);)

    auto dt = newTime - currentTime;

    for (auto& state : view)
    {
        auto const& [layout, ions, em] = state.ppc_update_populations();
        ionUpdater_.updatePopulations(ions, em, layout, dt, mode);
        // this needs to be done before calling the messenger
        view.model().resourcesManager->setTime(ions, *state.patch, newTime);
    }

    fromCoarser.fillIonGhostParticles(view.model().state.ions, level, newTime);
    fromCoarser.fillIonPopMomentGhosts(view.model().state.ions, level, newTime);

    for (auto& state : view)
    {
        ionUpdater_.updateIons(state.ions, state.layout);
        // no need to update time, since it has been done before
    }

    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    fromCoarser.fillIonMomentGhosts(view.model().state.ions, level, newTime);
}



} // namespace PHARE::solver




#endif
