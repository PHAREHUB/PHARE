
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

    // Hybrid
    using VecField_t   = typename core_types::VecField_t;
    using Grid_t       = typename core_types::Grid_t;
    using Electromag_t = typename core_types::Electromag_t;
    using Ions_t       = typename core_types::Ions_t;
    using GridLayout_t = typename core_types::GridLayout_t;
    using Electrons_t  = typename core_types::Electrons_t;

    // MHD
    using Grid_MHD       = typename core_types::Grid_MHD;
    using VecField_MHD   = typename core_types::VecField_MHD;
    using GridLayout_MHD = typename core_types::GridLayout_MHD;
    // core deps

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using HybridModel_t  = PHARE::solver::HybridModel<GridLayout_t, Electromag_t, Ions_t,
                                                      Electrons_t, PHARE::amr::SAMRAI_Types, Grid_t>;
    using MHDModel_t
        = PHARE::solver::MHDModel<GridLayout_MHD, VecField_MHD, PHARE::amr::SAMRAI_Types, Grid_MHD>;

    using SolverPPC_t = PHARE::solver::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;
    using SolverMHD_t
        = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types,
                                   typename decltype(opts)::template MHDTimeStepper_t<MHDModel_t>>;

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
