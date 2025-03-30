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
#include "core/numerics/ion_updater/new_updater.hpp"
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

    using ModelViews_t = HybridPPCModelView<HybridModel>;
    using Faraday_t    = typename ModelViews_t::Faraday_t;
    using Ampere_t     = typename ModelViews_t::Ampere_t;
    using Ohm_t        = typename ModelViews_t::Ohm_t;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

    Faraday_t faraday_;
    Ampere_t ampere_;
    Ohm_t ohm_;

    PHARE::core::other::IonUpdater<Ions, Electromag, GridLayout> ionUpdater_;


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



    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, ISolverModelView& views,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;


    void onRegrid() override { ionUpdater_.reset(); }


    std::shared_ptr<ISolverModelView> make_view(level_t& level, IPhysicalModel_t& model) override
    {
        return std::make_shared<ModelViews_t>(level, dynamic_cast<HybridModel&>(model));
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


    void saveState_(level_t& level, ModelViews_t& views);
    void restoreState_(level_t& level, ModelViews_t& views);


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


    // extend lifespan
    std::unordered_map<std::string, ParticleArray> tmpDomain;
    std::unordered_map<std::string, ParticleArray> patchGhost;

    template<typename Map>
    static void add_to(Map& map, std::string const& key, ParticleArray const& ps)
    {
        // vector copy drops the capacity (over allocation of the source)
        // we want to keep the overallocation somewhat - how much to be assessed
        ParticleArray empty{ps.box()};

        if (!map.count(key))
            map.emplace(key, empty);
        else
            map.at(key) = empty;

        auto& v = map.at(key);
        v.reserve(ps.capacity());
        v.replace_from(ps);
    }

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
void SolverPPC<HybridModel, AMR_Types>::saveState_(level_t& level, ModelViews_t& views)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::saveState_");

    for (auto& state : views)
    {
        std::stringstream ss;
        ss << state.patch->getGlobalId();
        for (auto& pop : state.ions)
        {
            std::string const key = ss.str() + "_" + pop.name();
            add_to(tmpDomain, key, pop.domainParticles());
            add_to(patchGhost, key, pop.patchGhostParticles());
        }
    }
}

template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::restoreState_(level_t& level, ModelViews_t& views)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::restoreState_");

    for (auto& state : views)
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
                                                     int const levelNumber, ISolverModelView& views,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::advanceLevel");

    auto& modelView   = dynamic_cast<ModelViews_t&>(views);
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
        fromCoarser.fillCurrentGhosts(views.model().state.J, level.getLevelNumber(), newTime);
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
        fromCoarser.fillCurrentGhosts(views.model().state.J, level.getLevelNumber(), newTime);
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
        fromCoarser.fillCurrentGhosts(views.model().state.J, levelNumber, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm");
        for (auto& state : views)
            state.electrons.update(state.layout);
        ohm_(views.layouts, views.N, views.Ve, views.Pe, views.electromag_B, views.J,
             views.electromag_E);
        setTime([](auto& state) -> auto& { return state.electromag.E; });

        fromCoarser.fillElectricGhosts(views.model().state.electromag.E, levelNumber, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, ModelViews_t& views,
                                                 Messenger& fromCoarser, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::average_");

    for (auto& state : views)
    {
        PHARE::core::average(state.electromag.B, state.electromagPred.B, state.electromagAvg.B);
        PHARE::core::average(state.electromag.E, state.electromagPred.E, state.electromagAvg.E);
    }

    // the following will fill E on all edges of all ghost cells, including those
    // on domain border. For level ghosts, electric field will be obtained from
    // next coarser level E average
    fromCoarser.fillElectricGhosts(electromagAvg_.E, level.getLevelNumber(), newTime);
}


template<typename... Args>
void _debug_log_move_ions(Args const&... args)
{
    auto const& [views] = std::forward_as_tuple(args...);

    std::size_t nbrDomainParticles        = 0;
    std::size_t nbrPatchGhostParticles    = 0;
    std::size_t nbrLevelGhostNewParticles = 0;
    std::size_t nbrLevelGhostOldParticles = 0;
    std::size_t nbrLevelGhostParticles    = 0; //
    for (auto& state : views)
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
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, ModelViews_t& views,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::moveIons_");
    PHARE_DEBUG_DO(_debug_log_move_ions(views);)

    TimeSetter setTime{views, newTime};

    {
        auto dt = newTime - currentTime;
        for (auto& state : views)
            ionUpdater_.updatePopulations(state.ions, state.electromagAvg, state.layout, dt, mode);
    }

    // this needs to be done before calling the messenger
    setTime([](auto& state) -> auto& { return state.ions; });

    fromCoarser.fillIonGhostParticles(views.model().state.ions, level, newTime);
    fromCoarser.fillIonPopMomentGhosts(views.model().state.ions, level, newTime);

    for (auto& state : views)
        ionUpdater_.updateIons(state.ions);
    // no need to update time, since it has been done before

    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    fromCoarser.fillIonMomentGhosts(views.model().state.ions, level, newTime);
}



} // namespace PHARE::solver




#endif
