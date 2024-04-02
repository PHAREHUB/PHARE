#ifndef PHARE_SOLVER_PIC_HPP
#define PHARE_SOLVER_PIC_HPP


#include "core/def/phare_mpi.hpp"

#include <SAMRAI/hier/Patch.h>


#include "amr/messengers/pic_messenger.hpp"
#include "amr/messengers/pic_messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"


#include "core/numerics/fermion_updater/fermion_updater.hpp"
#include "core/numerics/ampere/maxwell_ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

namespace PHARE::solver
{

template<typename PICModel, typename AMR_Types>
class SolverPIC : public ISolver<AMR_Types>
{
private:
    using Electromag       = typename PICModel::electromag_type;
    using Ions             = typename PICModel::ions_type;
    using ParticleArray    = typename Ions::particle_array_type;
    using VecFieldT        = typename PICModel::vecfield_type;
    using GridLayout       = typename PICModel::gridlayout_type;
    using ResourcesManager = typename PICModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using PICMessenger     = amr::PICMessenger<PICModel>;

    Electromag electromagPred_{"EMPred"};
    Electromag electromagAvg_{"EMAvg"};

    // are these ever called elsewhere? maybe delete
    static constexpr auto dimension    = PICModel::dimension;
    static constexpr auto interp_order = PICModel::gridlayout_type::interp_order;

    PHARE::core::Faraday<GridLayout> faraday_;
    PHARE::core::MaxwellAmpere<GridLayout> maxwellampere_;
    PHARE::core::FermionUpdater<Fermions, Electromag, VecField, GridLayout> fermionUpdater_;


public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

    SolverPIC()
        : ISolver<AMR_Types>{"PICSolver"}
    {
    }


    virtual ~SolverPIC() = default;

    virtual std::string modelName() const override { return PICModel::model_name; }


    virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    virtual void registerResources(IPhysicalModel_t& model) override;


    virtual void allocate(IPhysicalModel_t& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const override;
                        

    virtual void advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy, int const levelNumber,
                              IPhysicalModel_t& model, IMessenger& fromCoarserMessenger,
                              double const currentTime, double const newTime) override;

private:
    using Messenger = amr::PICMessenger<PICModel>;

    void clearVecField_(VecField& vecField, Indexes const&... ijk);

    void restartJ_(level_t& level, PICModel& model, Messenger& fromCoarser, double const currentTime);
                   
    void MAMF_(level_t& level, PICModel& model, Messenger& fromCoarser, double const currentTime, double const newTime);

    void moveFermions_(level_t& level, Fermions& fermions, PICModel& model, Electromag& electromag,
                       ResourcesManager& rm, Messenger& fromCoarser, double const currentTime,
                       double const newTime);

    void average_(level_t& level, PICModel& model, Messenger& fromCoarser);

    void set_(level_t& level, PICModel& model, Messenger& fromCoarser);

}; // end SolverPIC


template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::registerResources(IPhysicalModel_t& model)
{
    auto& pmodel = dynamic_cast<PICModel&>(model);
    pmodel.resourcesManager->registerResources(electromagPred_);
    pmodel.resourcesManager->registerResources(electromagAvg_);
}



template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::allocate(IPhysicalModel_t& model,
                                                 SAMRAI::hier::Patch& patch,
                                                 double const allocateTime) const
{
    auto& pmodel = dynamic_cast<PICModel&>(model);
    pmodel.resourcesManager->allocate(electromagPred_, patch, allocateTime);
    pmodel.resourcesManager->allocate(electromagAvg_, patch, allocateTime);
}


// TODO: this stays for now, but should probably be edited
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    auto const& Eavg  = electromagAvg_.E;
    auto const& Bpred = electromagPred_.B;

    hybridInfo.ghostElectric.emplace_back(core::VecFieldNames{Eavg});
    hybridInfo.initMagnetic.emplace_back(core::VecFieldNames{Bpred});
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                                                  int const levelNumber, IPhysicalModel_t& model,
                                                  IMessenger& fromCoarserMessenger,
                                                  double const currentTime, double const newTime)
{
    // Restart current density

    // Interpolate E, B at particle positions: included in boris' move()
    // Push particles
    // Project current onto the grid using Esirkepov scheme: included in fermion_updater

    // Solve Faraday, MaxwellAmpere

    // Average B in time and set Bavg to B


    PHARE_LOG_SCOPE("SolverPPC::advanceLevel");

    auto& PICModel         = dynamic_cast<PICModel&>(model);
    auto& PICState         = PICModel.state;
    auto& fromCoarser      = dynamic_cast<PICMessenger&>(fromCoarserMessenger);
    auto& resourcesManager = *PICModel.resourcesManager;
    auto level             = hierarchy->getPatchLevel(levelNumber);


    restartJ_(*level, PICModel, fromCoarser, currentTime);

    MAMF_(*level, PICModel, fromCoarser, currentTime, newTime);

    moveFermions_(*level, PICState.fermions, PICModel, electromagAvg_, resourcesManager, fromCoarser, currentTime,
                  newTime);

    average_(level_t& level, PICModel& model, Messenger& fromCoarser);

    set_(level_t& level, PICModel& model, Messenger& fromCoarser);
}


// used in restartJ_
template<typename VecField, typename... Indexes>
void SolverPIC<PICModel, AMR_Types>::clearVecField_(VecField& vecField, Indexes const&... ijk)
{
    auto& XField = vecField(Component::X);
    auto& YField = vecField(Component::Y);
    auto& ZField = vecField(Component::Z);

    auto& XField(ijk...) = 0.;
    auto& YField(ijk...) = 0.;
    auto& ZField(ijk...) = 0.;
}

// Restart current density (there's probably a better way to do this)
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::restartJ_(level_t& level, PICModel& model,
                                              Messenger& fromCoarser, double const currentTime)
{
    using LayoutHolder<GridLayout>::layout_;
    auto& PICState         = model.state;
    auto& resourcesManager = model.resourcesManager;

    auto& J = PICState.J;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, J);
        layout_->evalOnBox(J, [&](auto&... args) mutable { clearVecField_(J, args...); });
    }
    // TODO: check whether this is necessary
    // fromCoarser.fillCurrentGhosts(J, level.getLevelNumber(), currentTime);
}


// Solve Maxell-Faraday and Maxwell-Ampere
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::MAMF_(level_t& level, PICModel& model, Messenger& fromCoarser, 
                                           double const currentTime, double const newTime)
{
    auto& PICState         = model.state;
    auto& resourcesManager = model.resourcesManager;
    auto dt                = newTime - currentTime;
    auto& electromag       = PICState.electromag;

    {
        PHARE_LOG_SCOPE("SolverPPC::MAMF_.maxwellampere");

        auto& B     = electromag.B;
        auto& E     = electromag.E;
        auto& Epred = electromagPred_.E;
        auto& J     = PICState.J;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Epred, B, E, J);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, maxwellampere_);
            maxwellampere_(B, E, J, Epred, dt);

            resourcesManager->setTime(Epred, *patch, newTime);
        }
        fromCoarser.fillCurrentGhosts(E, level.getLevelNumber(), newTime);
    }

    {
        PHARE_LOG_SCOPE("SolverPPC::MAMF_.faraday");

        auto& Bpred = electromagPred_.B;
        auto& Epred = electromagPred_.E;
        auto& B     = electromag.B;

        for (auto& patch : level)
        {
            auto _      = resourcesManager->setOnPatch(*patch, Bpred, B, Epred);
            auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
            auto __     = core::SetLayout(&layout, faraday_);
            faraday_(B, Epred, Bpred, dt);
            resourcesManager->setTime(Bpred, *patch, newTime);
        }
    }
}

/* 
Electrons are to be treated as an ion population for the time being. The population should be dynamically 
created with the level and deleted with it.
This assumes interpolating flux and density is unnecessary (since bulk velocity and electron density 
are used only in Ohm's law. 
*/
template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::moveFermions_(level_t& level, Fermions& fermions, PICModel& model,
                                                  Electromag& electromag, ResourcesManager& rm,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::moveFermions_");

    PHARE_DEBUG_DO(std::size_t nbrDomainParticles = 0; std::size_t nbrPatchGhostParticles = 0;
                   std::size_t nbrLevelGhostNewParticles                                  = 0;
                   std::size_t nbrLevelGhostOldParticles                                  = 0;
                   std::size_t nbrLevelGhostParticles = 0; for (auto& patch
                                                                : level) {
                       auto _ = rm.setOnPatch(*patch, fermions); // TODO check if this is ok or if ions and electrons be treated separately

                       for (auto& particletype : fermions)
                       {
                        for (auto& pop : particletype)
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
                       }
                   })


    auto dt        = newTime - currentTime;
    auto& PICState = model.state;
    auto& J        = PICState.J;
    auto& ions     = fermions.ions;

    for (auto& patch : level)
    {
        auto _ = rm.setOnPatch(*patch, electromag, fermions, J);

        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        fermionUpdater_.updatePopulations(fermions, electromag, J, layout, dt);

        // this needs to be done before calling the messenger
        rm.setTime(fermions, *patch, newTime);
    }

    // Only ion ghosts are filled since electrons don't exist on the coarser level
    fromCoarser.fillIonGhostParticles(ions, level, newTime);
    fromCoarser.fillIonPopMomentGhosts(ions, level, newTime);

    // now Ni and Vi are calculated we can fill pure ghost nodes
    // these were not completed by the deposition of patch and levelghost particles
    fromCoarser.fillIonMomentGhosts(ions, level, newTime);
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::average_(level_t& level, PICModel& model,
                                                 Messenger& fromCoarser)
{
    PHARE_LOG_SCOPE("SolverPPC::average_");

    auto& PICState      = model.state;
    auto& resourcesManager = model.resourcesManager;

    auto& Bpred = electromagPred_.B;
    auto& Bavg  = electromagAvg_.B;
    auto& B     = PICState.electromag.B;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, Bpred, Bavg, B);
        PHARE::core::average(B, Bpred, Bavg);
    }
}

template<typename PICModel, typename AMR_Types>
void SolverPIC<PICModel, AMR_Types>::set_(level_t& level, PICModel& model,
                                                 Messenger& fromCoarser)
{
    PHARE_LOG_SCOPE("SolverPPC::set_");

    auto& PICState      = model.state;
    auto& resourcesManager = model.resourcesManager;

    auto& Bavg  = electromagAvg_.B;
    auto& B     = PICState.electromag.B;

    for (auto& patch : level)
    {
        auto _ = resourcesManager->setOnPatch(*patch, Bavg, B);
        B = Bavg;
    }
}

} // namespace PHARE::solver

#endif