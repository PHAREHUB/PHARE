#ifndef PHARE_HYBRID_LEVEL_INITIALIZER_HPP
#define PHARE_HYBRID_LEVEL_INITIALIZER_HPP

#include "mpi/mpi_utils.hpp"

#include "core/errors.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/physical_models/models.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/solvers/solver_hybrid_field_evolvers.hpp"
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

        using Ampere_t        = FieldEvolverDispatchers<HybridModel>::Ampere_t;
        using Ohm_t           = OhmLevelTransformer<HybridModel>;
        using ParticleArray_t = HybridModel::particle_array_type;
        using Interpolating_t
            = core::Interpolating<dimension, interp_order, /*atomic_interp*/ false>;


        core::OhmInfo ohm_info;


        inline bool isRootLevel(int const levelNumber) const { return levelNumber == 0; }

    public:
        using model_type = HybridModel;

        explicit HybridLevelInitializer(core::OhmInfo const& info)
            : ohm_info{info}
        {
        }

        explicit HybridLevelInitializer(PHARE::initializer::PHAREDict const& dict)
            : ohm_info{core::OhmInfo::FROM(dict["algo"]["ohm"])}
        {
        }

        virtual void initialize(hierarchy_ptr const& hierarchy, int levelNumber,
                                level_ptr const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) override
        {
            Interpolating_t interpolate_;

            auto& hybridModel = static_cast<HybridModel&>(model);
            auto& level       = amr_types::getLevel(*hierarchy, levelNumber);

            auto& hybMessenger = dynamic_cast<HybridMessenger&>(messenger);

            PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level");
            if (isRegridding)
            {
                PHARE_LOG_LINE_STR("regriding level " + std::to_string(levelNumber));
                PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::regrid");
                messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
            }
            else
            {
                if (isRootLevel(levelNumber))
                {
                    PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::root_init");
                    model.initialize(level);
                    messenger.fillRootGhosts(model, level, initDataTime);
                }
                else
                {
                    PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::initLevel");
                    messenger.initLevel(model, level, initDataTime);
                }
            }

            auto& ions = hybridModel.state.ions;
            auto& rm   = *hybridModel.resourcesManager;

            // overlapping patches were each filled over their whole domain,
            // keep particles only on the patch owning their cell
            auto const neighbors = amr::makeSameLevelNeighbors(*hierarchy, levelNumber);
            for (auto& patch : rm.enumerate(level, ions))
                ions.eraseDomainParticlesIn(
                    amr::makeNonOwnedBoxesFor<GridLayoutT::dimension>(*patch, neighbors));

            // now all particles are here, we must compute moments.

            auto accessor = amr::make_model_level_accessor(level, hybridModel, ions);
            auto const tiled
                = [](auto& view) -> auto& { return std::get<0>(view.args).chargeDensity(); };

            // per tile for tiled layouts: no tile reads another's moments
            {
                PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::domain_deposit");
                amr::level_exec(accessor, tiled, [&](auto& view, auto const&... tile_idx) {
                    auto& [ions_p] = view.args;
                    core::resetMoments(ions_p, tile_idx...);
                    core::depositParticles(ions_p, view.layout, interpolate_, core::DomainDeposit{},
                                           tile_idx...);
                });
            }


            // at this point flux and density is computed for all pops
            // but nodes on ghost box overlaps are not complete because they lack
            // contribution of neighbor particles.
            // The following two calls will += flux and density on these overlaps.
            {
                PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::border_sums");
                hybMessenger.fillFluxBorders(ions, level, initDataTime);
                hybMessenger.fillDensityBorders(ions, level, initDataTime);
            }

            // the only remaning incomplete nodes are those next to and on level ghost layers
            // we now complete them by depositing levelghost particles
            {
                PHARE_LOG_SCOPE(
                    1, "HybridLevelInitializer::initialize_level::level_ghost_deposit_and_moments");
                amr::level_exec(accessor, tiled, [&](auto& view, auto const&... tile_idx) {
                    auto& [ions_p] = view.args;
                    if (!isRootLevel(levelNumber))
                        core::depositParticles(ions_p, view.layout, interpolate_,
                                               core::LevelGhostDeposit{}, tile_idx...);

                    // now all nodes are complete, the total ion moments can safely be computed.
                    ions_p.computeChargeDensity(tile_idx...);
                    ions_p.computeBulkVelocity(tile_idx...);
                    ions_p.computeChargeFlux(tile_idx...);
                });
            }
            {
                PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::fillIonBorders");
                hybMessenger.fillIonBorders(ions, level, initDataTime);
            }

            // moments are not time interpolated: their level ghost nodes come from the
            // deposit of level ghost particles above. Only E and B ghosts are, from the
            // Bold_ and Eold_ copies saved by prepareStep() (see end of this function).


            // J is never messenger-communicated: compute it locally via Ampere on
            // physical+1 (exactly the layer Ohm reads) on EVERY level, init and
            // regrid paths included. B ghosts are filled on all paths above, so the
            // first ghost layer of J is the curl of those B ghosts.
            TimeSetter setTime{level, hybridModel, initDataTime};

            auto& B = hybridModel.state.electromag.B;
            auto& J = hybridModel.state.J;

            {
                PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::ampere");
                Ampere_t{level, hybridModel}(B, J);
                setTime(J);
            }

            // E only needs computing on a fresh root level: it is regridded on
            // regrid paths and initLevel has done it on new finer levels
            if (!isRegridding and isRootLevel(levelNumber))
            {
                PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level::root_ohm");
                auto& E         = hybridModel.state.electromag.E;
                auto& electrons = hybridModel.state.electrons;

                ElectronsLevelTransformer{level, hybridModel}(electrons);

                Ohm_t{ohm_info, level, hybridModel}(B, J, E, electrons);
                setTime(E);
                hybMessenger.fillElectricGhosts(E, level, 0.);
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
