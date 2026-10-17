#ifndef PHARE_HYBRID_LEVEL_INITIALIZER_HPP
#define PHARE_HYBRID_LEVEL_INITIALIZER_HPP

#include "core/errors.hpp"
#include "mpi/mpi_utils.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/solvers/solver_hybrid_field_evolvers.hpp"
#include "amr/physical_models/models.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/level_initializer/level_initializer.hpp"

#include "initializer/data_provider.hpp"

#include <exception>

namespace PHARE
{
namespace solver
{
    template<typename HybridModel>
    class HybridLevelInitializer : public LevelInitializer<typename HybridModel::amr_types>
    {
        using amr_types                    = HybridModel::amr_types;
        using level_ptr                    = amr_types::level_ptr;
        using hierarchy_ptr                = amr_types::hierarchy_ptr;
        using IPhysicalModelT              = IPhysicalModel<amr_types>;
        using HybridMessenger              = amr::HybridMessenger<HybridModel>;
        using GridLayoutT                  = HybridModel::gridlayout_type;
        static constexpr auto dimension    = GridLayoutT::dimension;
        static constexpr auto interp_order = GridLayoutT::options.interp_order;

        using Ampere_t = FieldEvolverDispatchers<HybridModel>::Ampere_t;
        using Ohm_t    = OhmLevelTransformer<HybridModel>;
        core::OhmInfo ohm_info;

        inline bool isRootLevel(int const levelNumber) const { return levelNumber == 0; }

    public:
        using model_type = HybridModel;

        explicit HybridLevelInitializer(PHARE::initializer::PHAREDict const& dict)
            : ohm_info{core::OhmInfo::FROM(dict["algo"]["ohm"])}
        {
        }

        virtual void initialize(hierarchy_ptr const& hierarchy, int levelNumber,
                                level_ptr const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) override
        {
            core::Interpolator<dimension, interp_order> interpolate_;
            auto& hybridModel = static_cast<HybridModel&>(model);
            auto& state       = hybridModel.state;
            auto& level       = amr_types::getLevel(*hierarchy, levelNumber);

            auto& hybMessenger = dynamic_cast<HybridMessenger&>(messenger);

            PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level");
            if (isRegridding)
            {
                PHARE_LOG_LINE_STR("regriding level " + std::to_string(levelNumber));
                PHARE_LOG_START(3, "hybridLevelInitializer::initialize : regriding block");
                messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
                PHARE_LOG_STOP(3, "hybridLevelInitializer::initialize : regriding block");
            }
            else
            {
                if (isRootLevel(levelNumber))
                {
                    PHARE_LOG_START(3, "hybridLevelInitializer::initialize : root level init");
                    model.initialize(level);
                    messenger.fillRootGhosts(model, level, initDataTime);
                    PHARE_LOG_STOP(3, "hybridLevelInitializer::initialize : root level init");
                }
                else
                {
                    PHARE_LOG_START(3, "hybridLevelInitializer::initialize : initlevel");
                    messenger.initLevel(model, level, initDataTime);
                    PHARE_LOG_STOP(3, "hybridLevelInitializer::initialize : initlevel");
                }
            }

            auto ions_accessor = amr::make_model_level_accessor(level, hybridModel, state.ions);

            // overlapping patches were each filled over their whole domain,
            // keep particles only on the patch owning their cell
            auto const neighbours = amr::makeSameLevelNeighbors(*hierarchy, levelNumber);
            amr::onLevel(ions_accessor, [&](auto& view) {
                auto& [ions] = view.args;
                ions.eraseDomainParticlesIn(view.nonOwnedBoxes(neighbours));
            });

            // now all particles are here, we must compute moments.
            amr::onLevel(ions_accessor, [&](auto& view) {
                auto& [ions] = view.args;
                core::resetMoments(ions);
                core::depositParticles(ions, view.layout, interpolate_, core::DomainDeposit{});
            });

            // at this point flux and density is computed for all pops
            // but nodes on ghost box overlaps are not complete because they lack
            // contribution of neighbor particles.
            // The following two calls will += flux and density on these overlaps.
            hybMessenger.fillFluxBorders(state.ions, level, initDataTime);
            hybMessenger.fillDensityBorders(state.ions, level, initDataTime);

            // the only remaning incomplete nodes are those next to and on level ghost layers
            // we now complete them by depositing levelghost particles
            amr::onLevel(ions_accessor, [&](auto& view) {
                auto& [ions] = view.args;
                if (!isRootLevel(levelNumber))
                    core::depositParticles(ions, view.layout, interpolate_,
                                           core::LevelGhostDeposit{});

                // now all nodes are complete, the total ion moments
                // can safely be computed.
                ions.computeChargeDensity();
                ions.computeBulkVelocity();
            });
            hybMessenger.fillIonBorders(state.ions, level, initDataTime);

            // on level i>0, this relies on 'prepareStep' having been called on when
            // level i-1 was initialized (at the end of this function)
            // it seems SAMRAI does not call timeInterpolate() at this point although
            // moments need time interpolation. It probably knows that
            // we are at a sync time across levels and that the time interpolation
            // is not needed. But is still seems to use the messenger temporaries like
            // NiOld etc. so prepareStep() must be called, see end of the function.
            // - TODO more better comment(s)
            // hybMessenger.fillIonMomentGhosts(ions, level, initDataTime);


            // J is never messenger-communicated: compute it locally via Ampere on
            // physical+1 (exactly the layer Ohm reads) on EVERY level, init and
            // regrid paths included. B ghosts are filled on all paths above, so the
            // first ghost layer of J is the curl of those B ghosts.
            TimeSetter setTime{level, hybridModel, initDataTime};

            Ampere_t{level, hybridModel}(state.electromag.B, state.J);
            setTime(state.J);

            // E only needs computing on a fresh root level: it is regridded on
            // regrid paths and initLevel has done it on new finer levels
            if (!isRegridding and isRootLevel(levelNumber))
            {
                amr::onLevel(amr::make_model_level_accessor(level, hybridModel, state.electrons),
                             [&](auto& view) {
                                 auto& [electrons] = view.args;
                                 electrons.update(view.layout);
                             });

                Ohm_t{ohm_info, level, hybridModel}(state.electromag.B, state.J, state.electromag.E,
                                                    state.electrons);
                setTime(state.electromag.E);
                hybMessenger.fillElectricGhosts(state.electromag.E, level, 0.);
            }


            // quantities have been computed on the level, like the moments
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
