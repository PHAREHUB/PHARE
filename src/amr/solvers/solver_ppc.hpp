#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/Patch.h>


#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "amr/solvers/solver.hpp"
// #include "amr/solvers/solver_ppc_model_view.hpp" // TORM

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


template<typename GridLayout, typename ResourcesManager, typename PatchLevel>
struct FaradayDispatcher
{
    FaradayDispatcher(ResourcesManager& resourcesManager, PatchLevel& level)
        : level_{level}
        , resourcesManager_{resourcesManager}
    {
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
    {
        auto Bs    = resourcesManager_.views(B);
        auto Es    = resourcesManager_.views(E);
        auto Bnews = resourcesManager_.views(Bnew);

        std::size_t const view_nbr = Bs.size();

        for (auto iview = 0; iview < view_nbr; ++iview)
        {
            auto& Bx = std::get<0>(Bs)[iview];
            auto& By = std::get<1>(Bs)[iview];
            auto& Bz = std::get<2>(Bs)[iview];
            // needs to make a vecfield out of these...
            // use Phil's TensorField base thing...
            //
            // core::Faraday<GridLayout> faraday_; // needs a layout...
            // faraday_(B, E, Bnew, dt);
        }
    }
    PatchLevel& level_;
    ResourcesManager& resourcesManager_;
};


template<typename GridLayout, typename ResourcesManager, typename PatchLevel>
struct AmpereDispatcher
{
    AmpereDispatcher(ResourcesManager& resourcesManager, PatchLevel& level)
        : level_{level}
        , resourcesManager_{resourcesManager}
    {
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField& J)
    {
        auto Bs = resourcesManager_.views(B);
        auto Js = resourcesManager_.views(J);

        std::size_t const view_nbr = Bs.size();

        for (auto iview = 0; iview < view_nbr; ++iview)
        {
            auto& Bx = std::get<0>(Bs)[iview];
            auto& By = std::get<1>(Bs)[iview];
            auto& Bz = std::get<2>(Bs)[iview];

            auto& Jx = std::get<0>(Js)[iview];
            auto& Jy = std::get<1>(Js)[iview];
            auto& Jz = std::get<2>(Js)[iview];
            // needs to make a vecfield out of these...
            // use Phil's TensorField base thing...
            //
            // core::Faraday<GridLayout> faraday_; // needs a layout...
            // faraday_(B, E, Bnew, dt);
        }
    }
    PatchLevel& level_;
    ResourcesManager& resourcesManager_;
};


template<typename GridLayout, typename ResourcesManager, typename PatchLevel>
struct OhmDispatcher
{
    OhmDispatcher(initializer::PHAREDict const& dict, ResourcesManager& resourcesManager,
                  PatchLevel& level)
        : /*ohm_{dict},*/
        level_{level}
        , resourcesManager_{resourcesManager}
    {
    }

    template<typename VecField, typename Field>
    void operator()(Field const& N, VecField const& Ve, Field const& Pe, VecField const& B,
                    VecField const& J, VecField& Enew)
    {
        auto Ns    = resourcesManager_.views(N);
        auto Ves   = resourcesManager_.views(Ve);
        auto Pes   = resourcesManager_.views(Pe);
        auto Bs    = resourcesManager_.views(B);
        auto Js    = resourcesManager_.views(J);
        auto Enews = resourcesManager_.views(Enew);

        std::size_t const view_nbr = Bs.size();

        for (auto iview = 0; iview < view_nbr; ++iview)
        {
            auto& Bx = std::get<0>(Bs)[iview];
            auto& By = std::get<1>(Bs)[iview];
            auto& Bz = std::get<2>(Bs)[iview];

            auto& Jx = std::get<0>(Js)[iview];
            auto& Jy = std::get<1>(Js)[iview];
            auto& Jz = std::get<2>(Js)[iview];
            // needs to make a vecfield out of these...
            // use Phil's TensorField base thing...
            //
            // core::Faraday<GridLayout> faraday_; // needs a layout...
            // faraday_(B, E, Bnew, dt);
        }
    }
    PatchLevel& level_;
    ResourcesManager& resourcesManager_;
};


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

    // TORM
    // using ModelViews_t = HybridPPCModelView<HybridModel>;
    // using Faraday_t    = typename ModelViews_t::Faraday_t;
    // using Ampere_t     = typename ModelViews_t::Ampere_t;
    // using Ohm_t        = typename ModelViews_t::Ohm_t;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

private:
    using FaradayDispatcher_t = FaradayDispatcher<GridLayout, ResourcesManager, level_t>;
    using AmpereDispatcher_t  = AmpereDispatcher<GridLayout, ResourcesManager, level_t>;
    using OhmDispatcher_t     = OhmDispatcher<GridLayout, ResourcesManager, level_t>;

    FaradayDispatcher_t faraday_;
    AmpereDispatcher_t ampere_;
    OhmDispatcher_t ohm_;

    PHARE::core::IonUpdater<Ions, Electromag, GridLayout> ionUpdater_;


public:
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



    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber,
                      IPhysicalModel<AMR_Types>& model, IMessenger& fromCoarserMessenger,
                      double const currentTime, double const newTime) override;


    void onRegrid() override { ionUpdater_.reset(); }


    // TORM
    // std::shared_ptr<ISolverModelView> make_view(level_t& level, IPhysicalModel_t& model) override
    // {
    //     return std::make_shared<ModelViews_t>(level, dynamic_cast<HybridModel&>(model));
    // }
    //

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


    void saveState_(level_t& level, HybridModel& model);
    void restoreState_(level_t& level, HybridModel& model);

    // struct TimeSetter
    // {
    //     template<typename QuantityAccessor>
    //     void operator()(QuantityAccessor accessor)
    //     {
    //         for (auto& state : views)
    //             views.model().resourcesManager->setTime(accessor(state), *state.patch, newTime);
    //     }
    //
    //     // ModelViews_t& views;
    //     double newTime;
    // };


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
void SolverPPC<HybridModel, AMR_Types>::saveState_(level_t& level, HybridModel& model)
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
void SolverPPC<HybridModel, AMR_Types>::restoreState_(level_t& level, HybridModel& model)
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
                                                     int const levelNumber,
                                                     IPhysicalModel<AMR_Types>& model,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::advanceLevel");

    auto& hybridModel = dynamic_cast<HybridModel&>(model);
    auto& fromCoarser = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto level        = hierarchy.getPatchLevel(levelNumber);

    predictor1_(*level, hybridModel, fromCoarser, currentTime, newTime);

    average_(*level, hybridModel, fromCoarser, newTime);

    saveState_(*level, hybridModel);

    moveIons_(*level, hybridModel, fromCoarser, currentTime, newTime,
              core::UpdaterMode::domain_only);

    predictor2_(*level, hybridModel, fromCoarser, currentTime, newTime);

    average_(*level, hybridModel, fromCoarser, newTime);

    restoreState_(*level, hybridModel);

    moveIons_(*level, hybridModel, fromCoarser, currentTime, newTime, core::UpdaterMode::all);

    corrector_(*level, hybridModel, fromCoarser, currentTime, newTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_");

    model.resources_manager().setTime(model.state, level, newTime);

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.faraday");
        auto dt     = newTime - currentTime;
        auto& B     = model.state.electromag.B;
        auto& E     = model.state.electromag.E;
        auto& Bpred = electromagPred_.B;

        faraday_(B, E, Bpred, dt);
        model.resources_manager().setTime(electromagPred_.B, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ampere");
        auto& B = electromagPred_.B;
        auto& J = model.state.J;
        ampere_(B, J);
        model.resources_manager().setTime(J, level, newTime);
        fromCoarser.fillCurrentGhosts(model.state.J, level.getLevelNumber(), newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ohm");
        // TODO
        // for (auto& state : views)
        //     state.electrons.update(state.layout);
        auto& N  = model.state.electrons.density();
        auto& Ve = model.state.electrons.velocity();
        auto& Pe = model.state.electrons.pressure();
        auto& B  = electromagPred_.B;
        auto& J  = model.state.J;
        auto& E  = electromagPred_.E;
        ohm_(N, Ve, Pe, B, J, E);
        model.resources_manager().setTime(E, level, newTime);
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_");

    model.resources_manager().setTime(model.state, level, newTime);

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.faraday");
        auto dt     = newTime - currentTime;
        auto& B     = model.state.electromag.B;
        auto& E     = electromagAvg_.E;
        auto& Bpred = electromagPred_.B;
        faraday_(B, E, Bpred, dt);
        model.resources_manager().setTime(electromagPred_.B, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ampere");
        auto& B = electromagPred_.B;
        auto& J = model.state.J;
        ampere_(B, J);
        model.resources_manager().setTime(J, level, newTime);
        fromCoarser.fillCurrentGhosts(model.state.J, level.getLevelNumber(), newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ohm");
        // TODO
        // for (auto& state : views)
        //     state.electrons.update(state.layout);
        auto& N  = model.state.electrons.density();
        auto& Ve = model.state.electrons.velocity();
        auto& Pe = model.state.electrons.pressure();
        auto& B  = electromagPred_.B;
        auto& J  = model.state.J;
        auto& E  = electromagPred_.E;
        ohm_(N, Ve, Pe, B, J, E);
        model.resources_manager().setTime(E, level, newTime);
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, HybridModel& model,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::corrector_");

    auto levelNumber = level.getLevelNumber();
    model.resources_manager().setTime(model.state, level, newTime);

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.faraday");
        auto dt = newTime - currentTime;
        auto& B = model.state.electromag.B;
        auto& E = electromagAvg_.E;

        faraday_(B, E, B, dt);
        model.resources_manager().setTime(B, level, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ampere");
        auto& B = model.state.electromag.B;
        auto& J = model.state.J;
        ampere_(B, J);
        model.resources_manager().setTime(J, level, newTime);
        fromCoarser.fillCurrentGhosts(model.state.J, level.getLevelNumber(), newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm");
        // TODO
        // for (auto& state : views)
        //     state.electrons.update(state.layout);
        auto& N  = model.state.electrons.density();
        auto& Ve = model.state.electrons.velocity();
        auto& Pe = model.state.electrons.pressure();
        auto& B  = model.state.electromag.B;
        auto& J  = model.state.J;
        auto& E  = model.state.electromag.E;
        ohm_(N, Ve, Pe, B, J, E);
        model.resources_manager().setTime(E, level, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, HybridModel& model,
                                                 Messenger& fromCoarser, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::average_");

    // TODO
    // for (auto& state : views)
    // {
    //     PHARE::core::average(state.electromag.B, state.electromagPred.B, state.electromagAvg.B);
    //     PHARE::core::average(state.electromag.E, state.electromagPred.E, state.electromagAvg.E);
    // }

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
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, HybridModel& model,
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
