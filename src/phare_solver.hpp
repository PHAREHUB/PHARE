
#ifndef PHARE_SOLVER_INCLUDE_HPP
#define PHARE_SOLVER_INCLUDE_HPP

#include "phare_amr.hpp" // IWYU pragma: keep

#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver_ppc.hpp"
#include "amr/multiphysics_integrator.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/messengers/messenger_factory.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/level_initializer/level_initializer_factory.hpp"
#include "python3/mhd_resolver.hpp"

namespace PHARE::solver
{
template<auto opts>
struct PHARE_Types
{
    auto static constexpr dimension     = opts.dimension;
    auto static constexpr interp_order  = opts.interp_order;
    auto static constexpr nbRefinedPart = opts.nbRefinedPart;

    // core deps
    using core_types = PHARE::core::PHARE_Types<opts>;

    struct Hybrid
    {
        using VecField_t   = core_types::Hybrid::VecField_t;
        using Grid_t       = core_types::Hybrid::Grid_t;
        using Electromag_t = core_types::Hybrid::Electromag_t;
        using Ions_t       = core_types::Hybrid::Ions_t;
        using GridLayout_t = core_types::Hybrid::GridLayout_t;
        using Electrons_t  = core_types::Hybrid::Electrons_t;
        using Model_t      = HybridModel< //
            GridLayout_t, Electromag_t, Ions_t, Electrons_t, amr::SAMRAI_Types, Grid_t>;
    };

    struct MHD
    {
        using Grid_t       = core_types::MHD::Grid_t;
        using VecField_t   = core_types::MHD::VecField_t;
        using GridLayout_t = core_types::MHD::GridLayout_t;
        using Model_t      = MHDModel<GridLayout_t, VecField_t, amr::SAMRAI_Types, Grid_t>;
    };

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using HybridModel_t  = Hybrid::Model_t;
    using MHDModel_t     = MHD::Model_t;

    using SolverPPC_t = PHARE::solver::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;
    using SolverMHD_t
        = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types,
                                   typename MHDResolver<opts, MHDModel_t>::MHDTimeStepper_t>;

    using LevelInitializerFactory_t
        = PHARE::solver::LevelInitializerFactory<HybridModel_t, MHDModel_t>;

    // amr deps
    using amr_types        = PHARE::amr::PHARE_Types<opts>;
    using RefinementParams = amr_types::RefinementParams;

    using MessengerFactory // = amr/solver bidirectional dependency
        = PHARE::amr::MessengerFactory<MHDModel_t, HybridModel_t, RefinementParams>;
    // amr deps

    using MultiPhysicsIntegrator_t
        = MultiPhysicsIntegrator<MessengerFactory, LevelInitializerFactory_t,
                                 PHARE::amr::SAMRAI_Types>;
};

} // namespace PHARE::solver

#endif // PHARE_SOLVER_INCLUDE_HPP
