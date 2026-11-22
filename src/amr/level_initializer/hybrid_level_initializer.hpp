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
#include "amr/physical_models/physical_model.hpp"
#include "amr/level_initializer/level_initializer.hpp"

#include "initializer/data_provider.hpp"

#include <exception>
#include <stdexcept>

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

        // one function per context: level 0 made from the initial condition or rebalanced, a
        // finer level made from its coarser level or regridded
        virtual void initialize(hierarchy_ptr const& hierarchy, int levelNumber,
                                level_ptr const& oldLevel, IPhysicalModelT& model,
                                amr::IMessenger<IPhysicalModelT>& messenger, double initDataTime,
                                bool isRegridding) override
        {
            PHARE_LOG_SCOPE(1, "HybridLevelInitializer::initialize_level");
            auto& hybridModel  = static_cast<HybridModel&>(model);
            auto& hybMessenger = dynamic_cast<HybridMessenger&>(messenger);

            if (isRootLevel(levelNumber))
            {
                if (isRegridding)
                    rebalanceL0_(hierarchy, oldLevel, hybridModel, hybMessenger, initDataTime);
                else
                    initL0_(hierarchy, hybridModel, hybMessenger, initDataTime);
            }
            else
            {
                if (isRegridding)
                    regridRefined_(hierarchy, levelNumber, oldLevel, hybridModel, hybMessenger,
                                   initDataTime);
                else
                    initRefined_(hierarchy, levelNumber, hybridModel, hybMessenger, initDataTime);
            }
        }

    private:
        // level 0 from the initial condition: particles loaded, B and E ghosts, the moments
        // deposited, J, then E from Ohm's law
        void initL0_(hierarchy_ptr const& hierarchy, HybridModel& model, HybridMessenger& messenger,
                     double const time)
        {
            PHARE_LOG_SCOPE(3, "HybridLevelInitializer::initL0_");
            auto& level = amr_types::getLevel(*hierarchy, 0);
            model.initialize(level);
            messenger.fillRootGhosts(model, level, time);
            keepOwnedParticles_(hierarchy, 0, level, model);
            depositMoments_(0, level, model, messenger, time);
            computeJ_(level, model, time);
            computeE_(level, model, messenger, time);
            messenger.prepareStep(model, level, time);
        }

        // level 0 onto new owners: no cell is new to it, everything is moved, the moments too
        // (messenger.regrid). They are not deposited again: their totals are the ones
        // synchronized from level 1, a deposit of level 0's particles would replace them
        void rebalanceL0_(hierarchy_ptr const& hierarchy, level_ptr const& oldLevel,
                          HybridModel& model, HybridMessenger& messenger, double const time)
        {
            PHARE_LOG_SCOPE(3, "HybridLevelInitializer::rebalanceL0_");
            if (!oldLevel)
                throw std::runtime_error("HybridLevelInitializer: level 0 regridded without "
                                         "its old level");
            auto& level = amr_types::getLevel(*hierarchy, 0);
            messenger.regrid(hierarchy, 0, oldLevel, model, time);
            keepOwnedParticles_(hierarchy, 0, level, model);
            computeJ_(level, model, time);
            messenger.prepareStep(model, level, time);
        }

        // a new finer level from its coarser level: B, E and the particles refined (E: by
        // initLevel), its level ghost particles split, the moments deposited, J
        void initRefined_(hierarchy_ptr const& hierarchy, int const levelNumber, HybridModel& model,
                          HybridMessenger& messenger, double const time)
        {
            PHARE_LOG_SCOPE(3, "HybridLevelInitializer::initRefined_");
            auto& level = amr_types::getLevel(*hierarchy, levelNumber);
            messenger.initLevel(hierarchy, model, level, time);
            finishRefined_(hierarchy, levelNumber, level, model, messenger, time);
        }

        // a finer level onto new boxes: its data moved where the old level was, refined from
        // the coarser level elsewhere, the moments deposited, J
        void regridRefined_(hierarchy_ptr const& hierarchy, int const levelNumber,
                            level_ptr const& oldLevel, HybridModel& model,
                            HybridMessenger& messenger, double const time)
        {
            PHARE_LOG_SCOPE(3, "HybridLevelInitializer::regridRefined_");
            auto& level = amr_types::getLevel(*hierarchy, levelNumber);
            messenger.regrid(hierarchy, levelNumber, oldLevel, model, time);
            finishRefined_(hierarchy, levelNumber, level, model, messenger, time);
        }

        // a finer level's particles and fields in place (initLevel or regrid): the particles
        // kept on their owning patch, the moments deposited, J, the old copies for the next step
        void finishRefined_(hierarchy_ptr const& hierarchy, int const levelNumber, auto& level,
                            HybridModel& model, HybridMessenger& messenger, double const time)
        {
            keepOwnedParticles_(hierarchy, levelNumber, level, model);
            depositMoments_(levelNumber, level, model, messenger, time);
            fillIonMomentGhosts_(hierarchy, level, model, messenger, time);
            computeJ_(level, model, time);
            messenger.prepareStep(model, level, time);
        }

        // the level ghost particles complete the domain nodes, border included, but not the
        // node beyond the level that Ohm reads: the total density and E-node flux ghosts are
        // refined from the coarser level in space and time, like E and B
        void fillIonMomentGhosts_(hierarchy_ptr const& hierarchy, auto& level, HybridModel& model,
                                  HybridMessenger& messenger, double const time)
        {
            auto& ions = model.state.ions;
            {
                // the time interpolation weight comes from the destination's time stamp
                TimeSetter setIonsTime{level, model, time};
                setIonsTime(ions);
            }
            messenger.fillIonMomentGhosts(hierarchy, ions, level, time);
        }

        // overlapping patches were each filled over their whole domain, keep particles only on
        // the patch owning their cell
        void keepOwnedParticles_(hierarchy_ptr const& hierarchy, int const levelNumber, auto& level,
                                 HybridModel& model)
        {
            auto& ions           = model.state.ions;
            auto const neighbors = amr::makeSameLevelNeighbors(*hierarchy, levelNumber);
            for (auto& patch : model.resourcesManager->enumerate(level, ions))
                ions.eraseDomainParticlesIn(
                    amr::makeNonOwnedBoxesFor<GridLayoutT::dimension>(*patch, neighbors));
        }

        // the moments from the particles: the domain ones, the level ghost ones (held by every
        // patch whose ghost box covers them, deposited by one of them), the border sums, then
        // the totals. Moments are not time interpolated: their level ghost nodes come from the
        // level ghost particles
        void depositMoments_(int const levelNumber, auto& level, HybridModel& model,
                             HybridMessenger& messenger, double const time)
        {
            core::Interpolator<dimension, interp_order> interpolate;
            auto& ions = model.state.ions;
            auto& rm   = *model.resourcesManager;
            for (auto& patch : rm.enumerate(level, ions))
            {
                auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
                core::resetMoments(ions);
                core::depositParticles(ions, layout, interpolate, core::DomainDeposit{});
            }
            if (!isRootLevel(levelNumber))
                messenger.depositLevelGhostParticlesOld(ions, level);

            // nodes on ghost box overlaps lack the neighbours' contributions: += them
            messenger.fillFluxBorders(ions, level, time);
            messenger.fillDensityBorders(ions, level, time);

            // every domain node complete, the totals. The ghost node beyond a refined level
            // that Ohm reads is not: fillIonMomentGhosts_
            for (auto& patch : rm.enumerate(level, ions))
            {
                ions.computeChargeDensity();
                ions.computeBulkVelocity();
                ions.computeChargeFlux();
            }
        }

        // J is never communicated: Ampere on physical+1 (the layer Ohm reads) on every level and
        // path, from the B ghosts each path filled
        void computeJ_(auto& level, HybridModel& model, double const time)
        {
            TimeSetter setTime{level, model, time};
            Ampere_t{level, model}(model.state.electromag.B, model.state.J);
            setTime(model.state.J);
        }

        // E from Ohm's law, on a level 0 made from the initial condition only: moved by a
        // rebalance or a regrid, refined by initLevel elsewhere
        void computeE_(auto& level, HybridModel& model, HybridMessenger& messenger,
                       double const time)
        {
            TimeSetter setTime{level, model, time};
            auto& E         = model.state.electromag.E;
            auto& electrons = model.state.electrons;
            for (auto& patch : model.resourcesManager->enumerate(level, electrons))
                electrons.update(amr::layoutFromPatch<GridLayoutT>(*patch));
            Ohm_t{ohm_info, level, model}(model.state.electromag.B, model.state.J, E, electrons);
            setTime(E);
            messenger.fillElectricGhosts(E, level, 0.);
        }
    };
} // namespace solver
} // namespace PHARE

#endif
