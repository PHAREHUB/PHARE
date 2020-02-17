

#ifndef PHARE_HYBRID_LEVEL_INITIALIZE_H
#define PHARE_HYBRID_LEVEL_INITIALIZE_H

#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/messenger.h"
#include "amr/resources_manager/amr_utils.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/moments/moments.h"
#include "level_initializer.h"
#include "solver/physical_models/hybrid_model.h"
#include "solver/physical_models/physical_model.h"
#include "core/data/ions/ions.h"

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
        using GridLayoutT                  = typename HybridModel::gridLayout_type;
        static constexpr auto dimension    = GridLayoutT::dimension;
        static constexpr auto interp_order = GridLayoutT::interp_order;


        inline bool isRootLevel(int levelNumber) const { return levelNumber == 0; }

    public:
        virtual void initialize(std::shared_ptr<hierarchy_t> const& hierarchy, int levelNumber,
                                std::shared_ptr<level_t> const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) const override
        {
            core::Interpolator<dimension, interp_order> interpolate_;
            auto& hybridModel = static_cast<HybridModel&>(model);
            auto& level       = amr_types::getLevel(*hierarchy, levelNumber);

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

                core::setNansOnGhosts<std::decay_t<decltype(ions)>, GridLayoutT>(ions, layout);
                ions.computeDensity();
                ions.computeBulkVelocity();
            }
        }
    };
} // namespace solver
} // namespace PHARE

#endif
