#ifndef PHARE_HYBRID_LEVEL_INITIALIZER_HPP
#define PHARE_HYBRID_LEVEL_INITIALIZER_HPP

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/moments/moments.hpp"
#include "amr/level_initializer/level_initializer.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "initializer/data_provider.hpp"


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
        explicit HybridLevelInitializer(PHARE::initializer::PHAREDict const& dict)
            : ohm_{dict["algo"]["ohm"]}
        {
        }
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
                PHARE_LOG_START("hybridLevelInitializer::initialize : root level init");
                model.initialize(level);
                messenger.fillRootGhosts(model, level, initDataTime);
                PHARE_LOG_STOP("hybridLevelInitializer::initialize : root level init");
            }

            else
            {
                if (isRegridding)
                {
                    std::cout << "reriding level " << levelNumber << "\n";
                    PHARE_LOG_START("hybridLevelInitializer::initialize : regriding block");
                    messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
                    PHARE_LOG_STOP("hybridLevelInitializer::initialize : regriding block");
                }
                else
                {
                    PHARE_LOG_START("hybridLevelInitializer::initialize : initlevel");
                    messenger.initLevel(model, level, initDataTime);
                    PHARE_LOG_STOP("hybridLevelInitializer::initialize : initlevel");
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

                ions.computeDensity();
                ions.computeBulkVelocity();
            }
            hybMessenger.prepareStep(hybridModel, level, initDataTime);
            hybMessenger.fillIonMomentGhosts(hybridModel.state.ions, level, initDataTime);


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
