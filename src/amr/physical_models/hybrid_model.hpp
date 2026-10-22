#ifndef PHARE_HYBRID_MODEL_HPP
#define PHARE_HYBRID_MODEL_HPP


#include "core/def.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/boundary/boundary_manager.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"

#include "initializer/data_provider.hpp"

#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/resources_manager.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include <map>
#include <limits>
#include <string>

namespace PHARE::solver
{
/**
 * @brief The HybridModel class is a concrete implementation of a IPhysicalModel. The class
 * holds a HybridState and a ResourcesManager.
 */
template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
class HybridModel : public IPhysicalModel<AMR_Types>
{
public:
    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::options.interp_order;

    using Interface              = IPhysicalModel<AMR_Types>;
    using amr_types              = AMR_Types;
    using electrons_t            = Electrons;
    using patch_t                = AMR_Types::patch_t;
    using level_t                = AMR_Types::level_t;
    using physical_quantity_type = core::HybridQuantity;
    using gridlayout_type        = GridLayoutT;
    using electromag_type        = Electromag;
    using vecfield_type          = Electromag::vecfield_type;
    using field_type             = vecfield_type::field_type;
    using grid_type              = Grid_t;
    using ions_type              = Ions;
    using particle_array_type    = Ions::particle_array_type;
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, grid_type>;
    using boundary_manager_type  = core::BoundaryManager<gridlayout_type, field_type>;
    using ParticleInitializerFactory
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;

    static constexpr std::string_view model_type_name = "HybridModel";
    static inline std::string const model_name{model_type_name};


    core::HybridState<Electromag, Ions, Electrons> state;
    std::shared_ptr<resources_manager_type> resourcesManager;
    std::shared_ptr<boundary_manager_type> boundaryManager;


    void initialize(level_t& level) override;


    /**
     * @brief allocate uses the ResourcesManager to allocate HybridState physical quantities on
     * the given Patch at the given allocateTime
     */
    virtual void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
    }




    /**
     * @brief fillMessengerInfo describes which variables of the model are to be initialized or
     * filled at ghost nodes.
     */
    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    NO_DISCARD auto setOnPatch(patch_t& patch)
    {
        return resourcesManager->setOnPatch(patch, *this);
    }


    // The electron pressure closure may save a state in beginStep (e.g. the
    // polytropic s = Pe/N^gamma) that its computePressure then reads on ghost nodes.
    // beginStep only sets it on domain nodes, since moments are not complete on all
    // ghost nodes (next to a level border); its ghosts are filled here, from
    // neighbor patches and, at level borders, interpolated from the coarser level's
    // state (saved when the coarser level started its current step, not time
    // interpolated). Must be called once per step and level, after every patch's
    // beginStep and before any computePressure. Nothing to do for a closure
    // without such a state. Ohm only reads Pe and Ve on domain nodes, so neither
    // needs ghosts.
    //
    // Needs a real refine operator (not nullptr): on level > 0 a patch's ghost
    // region can reach into territory only covered by the coarser level,
    // which needs genuine interpolation, not a same-resolution copy. Needs a
    // non-overwrite-interior fill pattern so the schedule only touches ghost
    // nodes, not the domain values of this patch.
    void fillElectronStepStateGhosts(int const levelNumber, double const fillTime)
    {
        auto const name = state.electrons.stepStateName();
        if (!name)
            return;

        if (!hierarchy_)
            throw std::runtime_error(
                "Error - HybridModel::fillElectronStepStateGhosts called before setHierarchy()");

        if (not electronGhostAlgoDeclared_)
        {
            auto&& [id] = resourcesManager->getIDsList(*name);
            stepStateGhostAlgo_.algo_->registerRefine(id, id, id, fieldRefineOp_,
                                                      nonOverwriteInteriorFieldFillPattern_);
            electronGhostAlgoDeclared_ = true;
        }

        // the refine operators only write fine nodes that are NaN (PHARE's marker for
        // "to be filled from the coarser level", see composite_field_refiner.hpp): without
        // this reset, ghost nodes at level borders keep the value of their first fill
        resetGhostsToNaN_(levelNumber, *name);

        stepStateGhostAlgo_.getOrCreateSchedule(hierarchy_, levelNumber).fillData(fillTime);
    }

    // A brand new level never gets Pe_'s domain cells set at all: unlike Ve_
    // (recomputed statelessly from ions/J every step, so it self-heals), Pe_
    // carries genuinely evolved state, like B/E. oldLevel is null when there
    // is no same-level source yet (initial refinement), or the replaced level
    // on regrid, whose overlap should be preserved rather than
    // re-interpolated from the coarser level.
    void initElectronPressureOnNewLevel(int const levelNumber, double const fillTime,
                                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel
                                        = nullptr)
    {
        if (!hierarchy_)
            throw std::runtime_error(
                "Error - HybridModel::initElectronPressureOnNewLevel called before setHierarchy()");

        auto&& [pe_id] = resourcesManager->getIDsList(state.electrons.pressure().name());

        SAMRAI::xfer::RefineAlgorithm initAlgo;
        initAlgo.registerRefine(pe_id, pe_id, pe_id, fieldRefineOp_);

        auto const level = hierarchy_->getPatchLevel(levelNumber);
        initAlgo.createSchedule(level, oldLevel, levelNumber - 1, hierarchy_)->fillData(fillTime);
    }

    // SolverPPC::advanceLevel only gets a bare hierarchy reference; only
    // HybridLevelInitializer::initialize() has the actual shared_ptr that
    // SAMRAI's coarse-fine createSchedule() needs, so it's cached here.
    void setHierarchy(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy)
    {
        hierarchy_ = hierarchy;
    }


    HybridModel(PHARE::initializer::PHAREDict const& dict,
                std::shared_ptr<resources_manager_type> const& _resourcesManager)
        : IPhysicalModel<AMR_Types>{model_name}
        , state{dict}
        , resourcesManager{_resourcesManager}
    {
        // no boundary conditions in hybrid for now, but boundaryManager needed by the messenger, so
        // we create a dummy one
        if (dict.contains("grid"))
            core::validatePhysicalBoundariesDeclared<dimension>(dict["grid"]);
        auto const has_bcs = dict.contains("grid") && dict["grid"].contains("domain_boundaries");
        boundaryManager    = std::make_shared<boundary_manager_type>(
            has_bcs ? dict["grid"]["domain_boundaries"] : PHARE::initializer::PHAREDict{},
            std::vector<core::HybridQuantity::Scalar>{},
            std::vector<core::HybridQuantity::Vector>{});
    }


    virtual ~HybridModel() override {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(state); }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;

private:
    // sets every non domain node of the named field to NaN, on every patch of the level
    void resetGhostsToNaN_(int const levelNumber, std::string const& name)
    {
        using FieldData_t = amr::FieldData<GridLayoutT, Grid_t>;

        auto&& [id] = resourcesManager->getIDsList(name);
        for (auto& patch : *hierarchy_->getPatchLevel(levelNumber))
        {
            auto& field       = FieldData_t::getField(*patch, id);
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);

            std::array<std::pair<std::uint32_t, std::uint32_t>, dimension> domain;
            for (std::size_t d = 0; d < dimension; ++d)
            {
                auto const [start, end] = layout.physicalStartToEnd(
                    field.physicalQuantity(), static_cast<core::Direction>(d));
                domain[d] = {start, end};
            }

            layout.evalOnGhostBox(field, [&](auto const&... ijk) {
                auto const index = std::array{ijk...};
                for (std::size_t d = 0; d < dimension; ++d)
                    if (index[d] < domain[d].first or index[d] > domain[d].second)
                    {
                        field(index) = std::numeric_limits<double>::quiet_NaN();
                        return;
                    }
            });
        }
    }

    // schedules are cached by level number, but a regrid replaces the
    // PatchLevel object for that number -- the weak_ptrs detect that so a
    // stale schedule isn't silently reused against a destroyed level. The
    // schedule also binds the next coarser level (level border ghosts are refined
    // from it), so a regrid of that level must invalidate it too.
    struct GhostFillAlgo
    {
        auto& getOrCreateSchedule(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                  int const ilvl)
        {
            auto const level           = hierarchy->getPatchLevel(ilvl);
            auto const coarser         = ilvl > 0 ? hierarchy->getPatchLevel(ilvl - 1) : nullptr;
            auto schedule_iter         = schedules_.find(ilvl);
            auto const create_schedule = schedule_iter == schedules_.end()
                                         or schedule_iter->second.level.lock() != level
                                         or schedule_iter->second.coarser.lock() != coarser;

            if (create_schedule)
                schedule_iter
                    = schedules_
                          .insert_or_assign(
                              ilvl, Entry{level, coarser,
                                          algo_->createSchedule(level, ilvl - 1, hierarchy)})
                          .first;

            return *schedule_iter->second.schedule;
        }

        struct Entry
        {
            std::weak_ptr<SAMRAI::hier::PatchLevel> level; // invalidated if the Level is destroyed
            std::weak_ptr<SAMRAI::hier::PatchLevel> coarser; // same for the next coarser level
            std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule;
        };

        std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> algo_
            = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
        std::map<int, Entry> schedules_;
    };

    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;

    bool electronGhostAlgoDeclared_ = false;
    GhostFillAlgo stepStateGhostAlgo_;

    // same kernel as the messenger's field refine operator; linear (order 2) is the only
    // refinement order supported (see amr::RefinementConfig)
    std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_
        = std::make_shared<amr::KernelFieldRefineOperator<GridLayoutT, Grid_t>>(
            amr::makeRefineKernel<GridLayoutT, Grid_t>(amr::FieldRefinementOrder::Linear));

    std::shared_ptr<amr::FieldFillPattern<dimension>> nonOverwriteInteriorFieldFillPattern_
        = std::make_shared<amr::FieldFillPattern<dimension>>();
};




//-------------------------------------------------------------------------
//                             definitions
//-------------------------------------------------------------------------


template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types, Grid_t>::initialize(
    level_t& level)
{
    auto& rm = *this->resourcesManager;
    for (auto& patch : rm.enumerate(level, state))
    {
        auto const layout = amr::layoutFromPatch<gridlayout_type>(*patch);

        for (auto& pop : state.ions)
            ParticleInitializerFactory::create(pop.particleInitializerInfo())
                ->loadParticles(pop.domainParticles(), layout);


        state.electrons.initialize(layout);
        state.electromag.initialize(layout);
    }
}



template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types, Grid_t>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    // only the charge density and the charge flux, which give the electron velocity,
    // are communicated. The mass density and the bulk velocity are diagnostics only
    hybridInfo.modelMagnetic      = state.electromag.B.name();
    hybridInfo.modelElectric      = state.electromag.E.name();
    hybridInfo.modelIonDensity    = state.ions.chargeDensityName();
    hybridInfo.modelIonChargeFlux = state.ions.chargeFlux().name();
    // electron pressure is a state for closures that evolve it (polytropic): the
    // coarser level must take the finer one's where they overlap, consistently
    // with the synchronized density, as for E
    hybridInfo.modelElectronPressure = state.electrons.pressure().name();

    hybridInfo.initElectric.emplace_back(state.electromag.E.name());
    hybridInfo.initMagnetic.emplace_back(state.electromag.B.name());

    hybridInfo.ghostElectric.push_back(hybridInfo.modelElectric);
    hybridInfo.ghostMagnetic.push_back(hybridInfo.modelMagnetic);

    auto transform_ = [](auto& ions, auto& inserter) {
        std::transform(std::begin(ions), std::end(ions), std::back_inserter(inserter),
                       [](auto const& pop) { return pop.name(); });
    };
    transform_(state.ions, hybridInfo.interiorParticles);
    transform_(state.ions, hybridInfo.levelGhostParticlesOld);
    transform_(state.ions, hybridInfo.levelGhostParticlesNew);
    transform_(state.ions, hybridInfo.patchGhostParticles);

    for (auto const& pop : state.ions)
    {
        hybridInfo.ghostFlux.emplace_back(pop.flux().name());
        hybridInfo.sumBorderFields.emplace_back(pop.particleDensity().name());
        hybridInfo.sumBorderFields.emplace_back(pop.chargeDensity().name());
    }

    hybridInfo.maxBorderFields.emplace_back(state.ions.chargeDensity().name());
    hybridInfo.maxBorderVecFields.emplace_back(state.ions.chargeFlux().name());
}




template<typename Model>
auto constexpr is_hybrid_model(Model* m) -> decltype(m->model_type_name, bool())
{
    return Model::model_type_name == "HybridModel";
}

template<typename... Args>
auto constexpr is_hybrid_model(Args...)
{
    return false;
}

template<typename Model>
auto constexpr is_hybrid_model_v = is_hybrid_model(static_cast<Model*>(nullptr));



} // namespace PHARE::solver

#endif
