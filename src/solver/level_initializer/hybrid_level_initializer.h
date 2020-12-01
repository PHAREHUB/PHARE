

#ifndef PHARE_HYBRID_LEVEL_INITIALIZER_H
#define PHARE_HYBRID_LEVEL_INITIALIZER_H

#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/messenger.h"
#include "amr/resources_manager/amr_utils.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/moments/moments.h"
#include "level_initializer.h"
#include "solver/physical_models/hybrid_model.h"
#include "solver/physical_models/physical_model.h"
#include "core/data/ions/ions.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/numerics/ohm/ohm.h"
#include "core/numerics/ampere/ampere.h"

namespace PHARE
{
namespace solver
{
    template<typename HybridModel>
    class HybridLevelInitializer : public LevelInitializer<typename HybridModel::amr_types>
    {
        using amr_types                    = typename HybridModel::amr_types;
        using hierarchy_t                  = typename amr_types::hierarchy_t;
        using level_t                      = typename amr_types::level_t;
        using patch_t                      = typename amr_types::patch_t;
        using IPhysicalModelT              = IPhysicalModel<amr_types>;
        using IMessengerT                  = amr::IMessenger<IPhysicalModelT>;
        using HybridMessenger              = amr::HybridMessenger<HybridModel>;
        using GridLayoutT                  = typename HybridModel::gridlayout_type;
        static constexpr auto dimension    = GridLayoutT::dimension;
        static constexpr auto interp_order = GridLayoutT::interp_order;


        PHARE::core::Ohm<GridLayoutT> ohm_;
        PHARE::core::Ampere<GridLayoutT> ampere_;


        inline bool isRootLevel(int levelNumber) const { return levelNumber == 0; }

    public:
        virtual void initialize(std::shared_ptr<hierarchy_t> const& hierarchy, int levelNumber,
                                std::shared_ptr<level_t> const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) override
        {
            core::Interpolator<dimension, interp_order> interpolate_;
            auto& hybridModel = static_cast<HybridModel&>(model);
            auto& level       = amr_types::getLevel(*hierarchy, levelNumber);

            auto& hybMessenger = dynamic_cast<HybridMessenger&>(messenger);


            if (isRootLevel(levelNumber))
            {
                model.initialize(level);
                messenger.fillRootGhosts(model, level, initDataTime);
            }

            else
            {
                if (isRegridding)
                {
                    messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
                }
                else
                {
                    messenger.initLevel(model, level, initDataTime);
                    messenger.prepareStep(model, level, initDataTime);
                }
            }

            // now all particles are here

            for (auto& patch : level)
            {
                auto& ions             = hybridModel.state.ions;
                auto& resourcesManager = hybridModel.resourcesManager;
                auto dataOnPatch       = resourcesManager->setOnPatch(*patch, ions);
                auto layout            = amr::layoutFromPatch<GridLayoutT>(*patch);


                core::resetMoments(ions);
                core::depositParticles(ions, layout, interpolate_, core::DomainDeposit{});
                core::depositParticles(ions, layout, interpolate_, core::PatchGhostDeposit{});


                if (!isRootLevel(levelNumber))
                {
                    core::depositParticles(ions, layout, interpolate_, core::LevelGhostDeposit{});
                }

                core::fixMomentGhosts(ions, layout);
                ions.computeDensity();
                ions.computeBulkVelocity();
            }


            if (isRootLevel(levelNumber))
            {
                auto& B = hybridModel.state.electromag.B;
                auto& J = hybridModel.state.J;

                for (auto& patch : level)
                {
                    auto _      = hybridModel.resourcesManager->setOnPatch(*patch, B, J);
                    auto layout = PHARE::amr::layoutFromPatch<GridLayoutT>(*patch);
                    auto __     = core::SetLayout(&layout, ampere_);
                    ampere_(B, J);

                    hybridModel.resourcesManager->setTime(J, *patch, 0.);
                }
                hybMessenger.fillCurrentGhosts(J, levelNumber, 0.);



                auto& electrons = hybridModel.state.electrons;
                auto& E         = hybridModel.state.electromag.E;

                for (auto& patch : level)
                {
                    auto layout = PHARE::amr::layoutFromPatch<GridLayoutT>(*patch);
                    auto _ = hybridModel.resourcesManager->setOnPatch(*patch, B, E, J, electrons);
                    electrons.update(layout);
                    auto& Ve = electrons.velocity();
                    auto& Ne = electrons.density();
                    auto& Pe = electrons.pressure();
                    auto __  = core::SetLayout(&layout, ohm_);
                    ohm_(Ne, Ve, Pe, B, J, E);
                    hybridModel.resourcesManager->setTime(E, *patch, 0.);
                }

                hybMessenger.fillElectricGhosts(E, levelNumber, 0.);
            }
        }
    };
} // namespace solver
} // namespace PHARE

#endif
