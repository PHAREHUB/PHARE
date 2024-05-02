#ifndef PHARE_HYBRID_LEVEL_INITIALIZER_HPP
#define PHARE_HYBRID_LEVEL_INITIALIZER_HPP

#include "amr/level_initializer/level_initializer.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/ions/ions.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/ohm/ohm.hpp"
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
                PHARE_LOG_START(3, "hybridLevelInitializer::initialize : root level init");
                model.initialize(level);
                messenger.fillRootGhosts(model, level, initDataTime);
                PHARE_LOG_STOP(3, "hybridLevelInitializer::initialize : root level init");
            }

            else
            {
                if (isRegridding)
                {
                    std::cout << "regriding level " << levelNumber << "\n";
                    PHARE_LOG_START(3, "hybridLevelInitializer::initialize : regriding block");
                    messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
                    PHARE_LOG_STOP(3, "hybridLevelInitializer::initialize : regriding block");
                }
                else
                {
                    PHARE_LOG_START(3, "hybridLevelInitializer::initialize : initlevel");
                    messenger.initLevel(model, level, initDataTime);
                    PHARE_LOG_STOP(3, "hybridLevelInitializer::initialize : initlevel");
                }
            }

            // now all particles are here
            // we must compute moments.

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

            // on level i>0, this relies on 'prepareStep' having been called on when
            // level i-1 was initialized (at the end of this function)
            // it seems SAMRAI does not call timeInterpolate() at this point although
            // both moments and J need time interpolation. It probably knows that
            // we are at a sync time across levels and that the time interpolation
            // is not needed. But is still seems to use the messenger tempoeraries like
            // NiOld etc. so prepareStep() must be called, see end of the function.
            hybMessenger.fillIonMomentGhosts(hybridModel.state.ions, level, initDataTime);


            // now moments are known everywhere, compute J and E
            // via Ampere and Ohm
            // this only needs to be done for the root level
            // since otherwise initLevel has done it already

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

            // quantities have been computed on the level,like the moments and J
            // that we later in the code need to get on level ghost nodes via
            // space and TIME interpolation. We thus need to save current values
            // in "old" messenger temporaries.
            // NOTE :  this may probably be skipped for finest level since, TBC at some point
            hybMessenger.prepareStep(hybridModel, level, initDataTime);
        }
    };
} // namespace solver
} // namespace PHARE

#endif
