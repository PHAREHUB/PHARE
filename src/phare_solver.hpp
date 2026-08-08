
#ifndef PHARE_SOLVER_INCLUDE_HPP
#define PHARE_SOLVER_INCLUDE_HPP

#include "phare_amr.hpp" // IWYU pragma: keep

#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver_ppc.hpp"
#include "amr/multiphysics_integrator.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/messengers/messenger_factory.hpp"
#include "amr/messengers/hybrid_hybrid_messenger_strategy.hpp"
#include "amr/messengers/mhd_hybrid_messenger_strategy.hpp"
#include "amr/messengers/mhd_messenger.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/data/particles/refine/splitter.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"
#include "amr/level_initializer/level_initializer_factory.hpp"
#include "amr/level_initializer/hybrid_level_initializer.hpp"
#include "amr/level_initializer/mhd_level_initializer.hpp"
#include "amr/solvers/mhd_resolver.hpp"
#include "amr/resources_manager/resources_manager.hpp"

#include <memory>

namespace PHARE::solver
{
// One specialization per enabled-combination: a disabled model's core types are never passed
// to ResourcesManager, so no resource-user specializations get instantiated/registered for a
// model this build doesn't want. When both models are enabled, they share one ResourcesManager
// instance/type (see comment on PHARE_Types::ResourcesManager_t below).
template<auto opts, typename CoreTypes, bool hasHybrid, bool hasMHD>
struct ResourcesManagerSelector;

template<auto opts, typename CoreTypes>
struct ResourcesManagerSelector<opts, CoreTypes, true, false>
{
    using type = amr::ResourcesManager<opts, typename CoreTypes::Hybrid>;
};

template<auto opts, typename CoreTypes>
struct ResourcesManagerSelector<opts, CoreTypes, false, true>
{
    using type = amr::ResourcesManager<opts, typename CoreTypes::MHD>;
};

template<auto opts, typename CoreTypes>
struct ResourcesManagerSelector<opts, CoreTypes, true, true>
{
    using type = amr::ResourcesManager<opts, typename CoreTypes::MHD, typename CoreTypes::Hybrid>;
};


// Bool-specialized holders: compile-time model enabling, to avoid expanding templates and
// emitting binary symbols for a model the opts value does not ask for. The `false`
// specialization is empty, so hybrid- (or mhd-) disabled opts never name the corresponding
// Model_t/Solver_t/etc. as a type.
// The `enabled` default is on the primary template because partial specializations may not
// carry default template arguments.
template<auto opts, typename CoreTypes, typename ResMan_t, bool enabled = has_hybrid_v<opts>>
struct HybridStack;

template<auto opts, typename CoreTypes, typename ResMan_t>
struct HybridStack<opts, CoreTypes, ResMan_t, false>
{
};

template<auto opts, typename CoreTypes, typename ResMan_t>
struct HybridStack<opts, CoreTypes, ResMan_t, true>
{
    struct Types
    {
        using GridLayout_t       = CoreTypes::Hybrid::GridLayout_t;
        using Electromag_t       = CoreTypes::Hybrid::Electromag_t;
        using Ions_t             = CoreTypes::Hybrid::Ions_t;
        using Electrons_t        = CoreTypes::Hybrid::Electrons_t;
        using Grid_t             = CoreTypes::Hybrid::Grid_t;
        using amr_types          = amr::SAMRAI_Types;
        using ResourcesManager_t = ResMan_t;
    };

    using GridLayout_t = CoreTypes::Hybrid::GridLayout_t;
    using Model_t      = HybridModel<Types>;
    using Solver_t     = PHARE::solver::SolverPPC<Model_t, PHARE::amr::SAMRAI_Types>;

    using Splitter_t = PHARE::amr::Splitter<PHARE::core::DimConst<opts.dimension>,
                                            PHARE::core::InterpConst<opts.interp_order>,
                                            PHARE::core::RefinedParticlesConst<opts.nbRefinedPart>>;
    using RefinementParams_t
        = PHARE::amr::RefinementParams<typename CoreTypes::Hybrid::ParticleArray_t, Splitter_t>;

    using LevelInitializer_t = HybridLevelInitializer<Model_t>;
};

template<auto opts, typename CoreTypes, typename ResMan_t, bool enabled = has_mhd_v<opts>>
struct MHDStack;

template<auto opts, typename CoreTypes, typename ResMan_t>
struct MHDStack<opts, CoreTypes, ResMan_t, false>
{
};

template<auto opts, typename CoreTypes, typename ResMan_t>
struct MHDStack<opts, CoreTypes, ResMan_t, true>
{
    struct Types
    {
        using GridLayout_t       = CoreTypes::MHD::GridLayout_t;
        using VecField_t         = CoreTypes::MHD::VecField_t;
        using Grid_t             = CoreTypes::MHD::Grid_t;
        using amr_types          = amr::SAMRAI_Types;
        using ResourcesManager_t = ResMan_t;
    };

    using GridLayout_t = CoreTypes::MHD::GridLayout_t;
    using Model_t       = MHDModel<Types>;
    using Solver_t = PHARE::solver::SolverMHD<Model_t, PHARE::amr::SAMRAI_Types,
                                              typename MHDResolver<opts, Model_t>::MHDTimeStepper_t>;

    using LevelInitializer_t = MHDLevelInitializer<Model_t>;
};

// One specialization per enabled-combination, so a disabled model's types are never named as
// template arguments to MessengerFactory/LevelInitializerFactory.
template<auto opts, typename Hybrid, typename MHD, bool hasHybrid, bool hasMHD>
struct FactorySelector;

// hybrid-only
template<auto opts, typename Hybrid, typename MHD>
struct FactorySelector<opts, Hybrid, MHD, true, false>
{
    using Messenger_t = amr::MessengerFactory<
        typename Hybrid::Model_t, typename Hybrid::Model_t,
        amr::HybridHybridMessengerStrategy<typename Hybrid::Model_t,
                                           typename Hybrid::RefinementParams_t>>;
    using LevelInit_t
        = LevelInitializerFactory<amr::SAMRAI_Types, typename Hybrid::LevelInitializer_t>;
};

// mhd-only
template<auto opts, typename Hybrid, typename MHD>
struct FactorySelector<opts, Hybrid, MHD, false, true>
{
    using Messenger_t = amr::MessengerFactory<typename MHD::Model_t, typename MHD::Model_t,
                                              amr::MHDMessenger<typename MHD::Model_t>>;
    using LevelInit_t = LevelInitializerFactory<amr::SAMRAI_Types, typename MHD::LevelInitializer_t>;
};

// both
template<auto opts, typename Hybrid, typename MHD>
struct FactorySelector<opts, Hybrid, MHD, true, true>
{
    using Messenger_t = amr::MessengerFactory<
        typename MHD::Model_t, typename Hybrid::Model_t,
        amr::HybridHybridMessengerStrategy<typename Hybrid::Model_t,
                                           typename Hybrid::RefinementParams_t>,
        amr::MHDHybridMessengerStrategy<typename MHD::Model_t, typename Hybrid::Model_t>,
        amr::MHDMessenger<typename MHD::Model_t>>;
    using LevelInit_t
        = LevelInitializerFactory<amr::SAMRAI_Types, typename Hybrid::LevelInitializer_t,
                                  typename MHD::LevelInitializer_t>;
};

template<auto opts>
struct PHARE_Types
{
    static_assert(opts.mhd_axes_consistent());
    static_assert(has_hybrid_v<opts> || has_mhd_v<opts>, "a build must enable at least one model");

    auto static constexpr dimension = opts.dimension;

    // core deps
    using core_types = PHARE::core::PHARE_Types<opts>;

    // shared by both models when both are enabled: one ResourcesManager instance can then serve
    // Hybrid and MHD resources - gated the same way as Hybrid/MHDStack below, so a disabled
    // model's core types are never passed to ResourcesManager (no resource-user specializations
    // instantiated/registered for a model this build doesn't want)
    using ResourcesManager_t =
        typename ResourcesManagerSelector<opts, core_types, has_hybrid_v<opts>,
                                          has_mhd_v<opts>>::type;

    using Hybrid = HybridStack<opts, core_types, ResourcesManager_t>;
    using MHD    = MHDStack<opts, core_types, ResourcesManager_t>;

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;

    // amr deps
    using amr_types = PHARE::amr::PHARE_Types<opts>;

    using Selector_t = FactorySelector<opts, Hybrid, MHD, has_hybrid_v<opts>, has_mhd_v<opts>>;

    using MessengerFactory // = amr/solver bidirectional dependency
        = Selector_t::Messenger_t;
    using LevelInitializerFactory_t = Selector_t::LevelInit_t;
    // amr deps

    using MultiPhysicsIntegrator_t
        = MultiPhysicsIntegrator<MessengerFactory, LevelInitializerFactory_t,
                                 PHARE::amr::SAMRAI_Types>;
};

// Same bool-specialization trick as HybridStack/MHDStack above, for runtime storage rather than
// types: the `false` specialization is empty, so a disabled model's resource manager / model
// shared_ptrs are simply not declared. Simulator never names the disabled model as a type outside
// a lazily-instantiated (non-virtual) function body.
template<auto opts, bool enabled = has_hybrid_v<opts>>
struct HybridSimState
{
};

template<auto opts>
struct HybridSimState<opts, true>
{
    using Model_t  = PHARE_Types<opts>::Hybrid::Model_t;
    using ResMan_t = Model_t::resources_manager_type;
    std::shared_ptr<ResMan_t> resman_;
    std::shared_ptr<Model_t> model_;
};

template<auto opts, bool enabled = has_mhd_v<opts>>
struct MHDSimState
{
};

template<auto opts>
struct MHDSimState<opts, true>
{
    using Model_t  = PHARE_Types<opts>::MHD::Model_t;
    using ResMan_t = Model_t::resources_manager_type;
    std::shared_ptr<ResMan_t> resman_;
    std::shared_ptr<Model_t> model_;
};

} // namespace PHARE::solver

#endif // PHARE_SOLVER_INCLUDE_HPP
