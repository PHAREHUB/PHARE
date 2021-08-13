#ifndef PHARE_SOLVER_PPC_GPU_H
#define PHARE_SOLVER_PPC_GPU_H


#include <cassert>

#include "amr/solvers/solver_ppc.h"
#include "core/numerics/ion_updater/gpu/ion_updater.h"

#include "amr/solvers/gpu/offloader.h"

namespace PHARE::solver::gpu
{
// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
class SolverPPC : public PHARE::solver::SolverPPC<HybridModel, AMR_Types>
{
public:
    static constexpr auto dimension    = HybridModel::dimension;
    static constexpr auto interp_order = HybridModel::gridlayout_type::interp_order;

    using Super            = PHARE::solver::SolverPPC<HybridModel, AMR_Types>;
    using This             = PHARE::solver::gpu::SolverPPC<HybridModel, AMR_Types>;
    using HybridModel_t    = HybridModel;
    using Electromag       = typename HybridModel::electromag_type;
    using Ions             = typename HybridModel::ions_type;
    using ParticleArray    = typename Ions::particle_array_type;
    using VecFieldT        = typename HybridModel::vecfield_type;
    using GridLayout       = typename HybridModel::gridlayout_type;
    using ResourcesManager = typename HybridModel::resources_manager_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using HybridMessenger  = amr::HybridMessenger<HybridModel>;

    using Offloader_t = PHARE::solver::gpu::Offloader<This>;

private:
    using Super::average_;
    using Super::corrector_;
    using Super::electromagAvg_;
    using Super::predictor1_;
    using Super::predictor2_;
    using Super::restoreState_;
    using Super::saveState_;

    //     PHARE::core::Faraday<GridLayout> faraday_;
    //     PHARE::core::Ampere<GridLayout> ampere_;
    //     PHARE::core::Ohm<GridLayout> ohm_;

    std::unique_ptr<Offloader_t> offloader;
    PHARE::core::gpu::IonUpdater<Ions, Electromag, GridLayout, Offloader_t> ionUpdater_;
    //     PHARE::core::IonUpdater<Ions, Electromag, GridLayout> ionUpdater_;


public:
    using patch_t     = typename AMR_Types::patch_t;
    using level_t     = typename AMR_Types::level_t;
    using hierarchy_t = typename AMR_Types::hierarchy_t;

    explicit SolverPPC(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
        , offloader{std::make_unique<PHARE::solver::gpu_mkn::Offloader<This>>(dict)}
        , ionUpdater_{dict["ion_updater"], *offloader}
    //     , ionUpdater_{dict["ion_updater"]}

    {
    }

    virtual ~SolverPPC() = default;


    virtual std::string modelName() const override { return HybridModel::model_name; }


    virtual void advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy, int const levelNumber,
                              IPhysicalModel_t& model, IMessenger& fromCoarserMessenger,
                              double const currentTime, double const newTime) override;

private:
    using Messenger = amr::HybridMessenger<HybridModel>;

    void moveIons_(level_t& level, Ions& ions, Electromag& electromag, ResourcesManager& rm,
                   Messenger& fromCoarser, double const currentTime, double const newTime,
                   core::UpdaterMode mode);

}; // end solverPPC


// -----------------------------------------------------------------------------

template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::advanceLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                                                     int const levelNumber, IPhysicalModel_t& model,
                                                     IMessenger& fromCoarserMessenger,
                                                     double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE("SolverPPC::advanceLevel");

    auto& hybridModel      = dynamic_cast<HybridModel&>(model);
    auto& hybridState      = hybridModel.state;
    auto& fromCoarser      = dynamic_cast<HybridMessenger&>(fromCoarserMessenger);
    auto& resourcesManager = *hybridModel.resourcesManager;
    auto level             = hierarchy->getPatchLevel(levelNumber);

    auto reset_moments = [&]() {
        auto& ions = hybridState.ions;
        for (auto& patch : *level)
        {
            auto _ = resourcesManager.setOnPatch(*patch, ions);
            resetMoments(ions);
        }
    };
    reset_moments();

    // GPU
    offloader->clear();
    for (auto& patch : *level)
    {
        auto _      = resourcesManager.setOnPatch(*patch, hybridState);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        offloader->alloc(layout, hybridState, newTime - currentTime);
    }
    // GPU

    predictor1_(*level, hybridModel, fromCoarser, currentTime, newTime);
    average_(*level, hybridModel);
    //     saveState_(*level, hybridState.ions, resourcesManager);

    // GPU
    offloader->template move<0>().join();
    // GPU

    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::domain_only);

    predictor2_(*level, hybridModel, fromCoarser, currentTime, newTime);
    average_(*level, hybridModel);
    //     restoreState_(*level, hybridState.ions, resourcesManager);

    reset_moments();

    // GPU
    offloader->template move<1>().join();
    // GPU

    moveIons_(*level, hybridState.ions, electromagAvg_, resourcesManager, fromCoarser, currentTime,
              newTime, core::UpdaterMode::all);

    corrector_(*level, hybridModel, fromCoarser, currentTime, newTime);
}




template<typename HybridModel, typename AMR_Types>
void SolverPPC<HybridModel, AMR_Types>::moveIons_(level_t& level, Ions& ions,
                                                  Electromag& electromag, ResourcesManager& rm,
                                                  Messenger& fromCoarser, double const currentTime,
                                                  double const newTime, core::UpdaterMode mode)
{
    PHARE_LOG_SCOPE("SolverPPC::moveIons_");

    std::size_t nbrDomainParticles        = 0;
    std::size_t nbrPatchGhostParticles    = 0;
    std::size_t nbrLevelGhostNewParticles = 0;
    std::size_t nbrLevelGhostOldParticles = 0;
    std::size_t nbrLevelGhostParticles    = 0;
    for (auto& patch : level)
    {
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
                throw std::runtime_error("Error - there are less old level ghost particles ("
                                         + std::to_string(nbrLevelGhostOldParticles)
                                         + ") than pushable ("
                                         + std::to_string(nbrLevelGhostParticles) + ")");
        }
    }


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
    fromCoarser.fillIonMomentGhosts(ions, level, currentTime, newTime);

    for (auto& patch : level)
    {
        auto _      = rm.setOnPatch(*patch, electromag, ions);
        auto layout = PHARE::amr::layoutFromPatch<GridLayout>(*patch);
        ionUpdater_.updateIons(ions, layout);

        // no need to update time, since it has been done before
    }
}
} // namespace PHARE::solver::gpu


#endif
