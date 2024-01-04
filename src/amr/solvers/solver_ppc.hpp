#ifndef PHARE_SOLVER_PPC_HPP
#define PHARE_SOLVER_PPC_HPP

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/Patch.h>


#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "amr/solvers/solver.hpp"


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


    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};


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

    virtual ~SolverPPC() = default;


    virtual std::string modelName() const override { return HybridModel::model_name; }

    virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    virtual void registerResources(IPhysicalModel_t& model) override;


    virtual void allocate(IPhysicalModel_t& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override;



    virtual void advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy, int const levelNumber,
                              IPhysicalModel_t& model, IMessenger& fromCoarserMessenger,
                              double const currentTime, double const newTime) override;


    void onRegrid() override { ionUpdater_.reset(); }

private:
    using Messenger = amr::HybridMessenger<HybridModel>;


    void predictor1_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void predictor2_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                     double const currentTime, double const newTime);


    void corrector_(level_t& level, HybridModel& model, Messenger& fromCoarser,
                    double const currentTime, double const newTime);


    void average_(level_t& level, HybridModel& model, Messenger& fromCoarser, double const newTime);


    void moveIons_(level_t& level, Ions& ions, Electromag& electromag, ResourcesManager& rm,
                   Messenger& fromCoarser, double const currentTime, double const newTime,
                   core::UpdaterMode mode);


    void saveState_(level_t& level, Ions& ions, ResourcesManager& rm);

    void restoreState_(level_t& level, Ions& ions, ResourcesManager& rm);


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
void SolverPPC<HybridModel, AMR_Types>::saveState_(level_t& level, Ions& ions, ResourcesManager& rm)
{
    for (auto& patch : level)
    {
        std::stringstream ss;
        ss << patch->getGlobalId();

        auto _ = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
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
void SolverPPC<HybridModel, AMR_Types>::restoreState_(level_t& level, Ions& ions,
                                                      ResourcesManager& rm)
{
    for (auto& patch : level)
    {
        std::stringstream ss;
        ss << patch->getGlobalId();

        auto _ = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            pop.domainParticles()     = std::move(tmpDomain.at(ss.str() + "_" + pop.name()));
            pop.patchGhostParticles() = std::move(patchGhost.at(ss.str() + "_" + pop.name()));
        }
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                                                     int const levelNumber, IPhysicalModel_t& model,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::advanceLevel");

    auto& hybridModel      = dynamic_cast<HybridModel&>(model);
    auto& hybridState      = hybridModel.state;
    auto& fromCoarser      = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto& resourcesManager = *hybridModel.resourcesManager;
    auto level             = hierarchy->getPatchLevel(levelNumber);


    predictor1_(*level, hybridModel, fromCoarser, currentTime, newTime);

    average_(*level, hybridModel, fromCoarser, newTime);

    saveState_(*level, hybridState.ions, resourcesManager);

    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::domain_only);

    predictor2_(*level, hybridModel, fromCoarser, currentTime, newTime);


    average_(*level, hybridModel, fromCoarser, newTime);

    restoreState_(*level, hybridState.ions, resourcesManager);

    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::all);

    corrector_(*level, hybridModel, fromCoarser, currentTime, newTime);


    // return newTime;
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor1_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto& electromag       = hybridState.electromag;
    auto levelNumber       = level.getLevelNumber();


    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.faraday");

        auto& Bpred = electromagPred_.B;
        auto& B     = electromag.B;
        auto& E     = electromag.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, B, E);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, E, Bpred, dt);
            resourcesManager->setTime(Bpred, *patch, newTime);
        }
    }



    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ampere");

        auto& Bpred = electromagPred_.B;
        auto& J     = hybridState.J;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, ampere_);
            ampere_(Bpred, J);

            resourcesManager->setTime(J, *patch, newTime);
        }
        fromCoarser.fillCurrentGhosts(J, levelNumber, newTime);
    }



    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor1_.ohm");

        auto& electrons = hybridState.electrons;
        auto& Bpred     = electromagPred_.B;
        auto& Epred     = electromagPred_.E;
        auto& J         = hybridState.J;

        for (auto& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, Epred, J, electrons);
            electrons.update(layout);
            auto& Ve = electrons.velocity();
            auto& Ne = electrons.density();
            auto& Pe = electrons.pressure();
            auto __  = core::SetLayout(&layout, ohm_);
            ohm_(Ne, Ve, Pe, Bpred, J, Epred);
            resourcesManager->setTime(Epred, *patch, newTime);
        }
    }
}


template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::predictor2_(level_t& level, HybridModel& model,
                                                    Messenger& fromCoarser,
                                                    double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto levelNumber       = level.getLevelNumber();



    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.faraday");

        auto& Bpred = electromagPred_.B;
        auto& B     = hybridState.electromag.B;
        auto& Eavg  = electromagAvg_.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, B, Eavg);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, Eavg, Bpred, dt);

            resourcesManager->setTime(Bpred, *patch, newTime);
        }
    }


    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ampere");

        auto& Bpred = electromagPred_.B;
        auto& J     = hybridState.J;


        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, ampere_);
            ampere_(Bpred, J);

            resourcesManager->setTime(J, *patch, newTime);
        }
        fromCoarser.fillCurrentGhosts(J, levelNumber, newTime);
    }


    {
        PHARE_LOG_SCOPE(1, "SolverPPC::predictor2_.ohm");

        auto& electrons = hybridState.electrons;
        auto& Bpred     = electromagPred_.B;
        auto& Epred     = electromagPred_.E;
        auto& J         = hybridState.J;

        for (auto& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, Epred, J, electrons);
            electrons.update(layout);
            auto& Ve = electrons.velocity();
            auto& Ne = electrons.density();
            auto& Pe = electrons.pressure();
            auto __  = core::SetLayout(&layout, ohm_);
            ohm_(Ne, Ve, Pe, Bpred, J, Epred);
            resourcesManager->setTime(Epred, *patch, newTime);
        }
    }
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::corrector_(level_t& level, HybridModel& model,
                                                   Messenger& fromCoarser, double const currentTime,
                                                   double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::corrector_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto levelNumber       = level.getLevelNumber();

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.faraday");

        auto& B    = hybridState.electromag.B;
        auto& Eavg = electromagAvg_.E;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, B, Eavg);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, Eavg, B, dt);

            resourcesManager->setTime(B, *patch, newTime);
        }
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ampere");

        auto& B = hybridState.electromag.B;
        auto& J = hybridState.J;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, B, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, ampere_);
            ampere_(B, J);

            resourcesManager->setTime(J, *patch, newTime);
        }
        fromCoarser.fillCurrentGhosts(J, levelNumber, newTime);
    }

    {
        PHARE_LOG_SCOPE(1, "SolverPPC::corrector_.ohm");

        auto& electrons = hybridState.electrons;
        auto& B         = hybridState.electromag.B;
        auto& E         = hybridState.electromag.E;
        auto& J         = hybridState.J;

        for (auto& patch : level)
        {
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, B, E, J, electrons);
            electrons.update(layout);
            auto& Ve = electrons.velocity();
            auto& Ne = electrons.density();
            auto& Pe = electrons.pressure();
            auto __  = core::SetLayout(&layout, ohm_);
            ohm_(Ne, Ve, Pe, B, J, E);
            resourcesManager->setTime(E, *patch, newTime);
        }

        fromCoarser.fillElectricGhosts(E, levelNumber, newTime);
    }
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::average_(level_t& level, HybridModel& model,
                                                 Messenger& fromCoarser, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::average_");

    auto& hybridState      = model.state;
    auto& resourcesManager = model.resourcesManager;

    auto& Epred = electromagPred_.E;
    auto& Bpred = electromagPred_.B;
    auto& Bavg  = electromagAvg_.B;
    auto& Eavg  = electromagAvg_.E;
    auto& B     = hybridState.electromag.B;
    auto& E     = hybridState.electromag.E;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, electromagAvg_, electromagPred_,
                                              hybridState.electromag);
        PHARE::core::average(B, Bpred, Bavg);
        PHARE::core::average(E, Epred, Eavg);
    }
    // the following will fill E on all edges of all ghost cells, including those
    // on domain border. For level ghosts, electric field will be obtained from
    // next coarser level E average
    fromCoarser.fillElectricGhosts(Eavg, level.getLevelNumber(), newTime);
}



template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, Ions& ions,
                                                  Electromag& electromag, ResourcesManager& rm,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    PHARE_LOG_SCOPE(1, "SolverPPC::moveIons_");

    PHARE_DEBUG_DO(std::size_t nbrDomainParticles = 0; std::size_t nbrPatchGhostParticles = 0;
                   std::size_t nbrLevelGhostNewParticles                                  = 0;
                   std::size_t nbrLevelGhostOldParticles                                  = 0;
                   std::size_t nbrLevelGhostParticles = 0; for (auto& patch
                                                                : level) {
                       auto _ = rm.setOnPatch(*patch, ions);

                       for (auto& pop : ions)
                       {
                           nbrDomainParticles += pop.domainParticles().size();
                           nbrPatchGhostParticles += pop.patchGhostParticles().size();
                           nbrLevelGhostNewParticles += pop.levelGhostParticlesNew().size();
                           nbrLevelGhostOldParticles += pop.levelGhostParticlesOld().size();
                           nbrLevelGhostParticles += pop.levelGhostParticles().size();
                           nbrPatchGhostParticles += pop.patchGhostParticles().size();

                           if (nbrLevelGhostOldParticles < nbrLevelGhostParticles
                               and nbrLevelGhostOldParticles > 0)
                               throw std::runtime_error(
                                   "Error - there are less old level ghost particles ("
                                   + std::to_string(nbrLevelGhostOldParticles) + ") than pushable ("
                                   + std::to_string(nbrLevelGhostParticles) + ")");
                       }
                   })


    auto dt = newTime - currentTime;

    for (auto& patch : level)
    {
        auto _ = rm.setOnPatch(*patch, electromag, ions);

        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        ionUpdater_.updatePopulations(ions, electromag, layout, dt, mode);

        // this needs to be done before calling the messenger
        rm.setTime(ions, *patch, newTime);
    }

    fromCoarser.fillIonGhostParticles(ions, level, newTime);
    fromCoarser.fillIonPopMomentGhosts(ions, level, newTime);

    for (auto& patch : level)
    {
        auto _      = rm.setOnPatch(*patch, electromag, ions);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        ionUpdater_.updateIons(ions, layout);

        // no need to update time, since it has been done before
    }
    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    fromCoarser.fillIonMomentGhosts(ions, level, newTime);
}


} // namespace PHARE::solver


#endif
