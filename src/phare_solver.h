
#ifndef PHARE_SOLVER_INCLUDE_H
#define PHARE_SOLVER_INCLUDE_H

#include "phare_amr.h"

#include "amr/solvers/solver.h"
#include "amr/solvers/solver_mhd.h"
#include "amr/solvers/solver_ppc.h"

#include "amr/level_initializer/level_initializer.h"
#include "amr/level_initializer/level_initializer_factory.h"
#include "amr/multiphysics_integrator.h"
#include "amr/physical_models/hybrid_model.h"
#include "amr/physical_models/mhd_model.h"
#include "amr/physical_models/physical_model.h"

#if defined(PHARE_WITH_GPU)
#include "amr/solvers/gpu/solver_ppc.h"
#elif defined(WITH_RAJA) and defined(WITH_UMPIRE)
#include "amr/solvers/llnl/solver_ppc.h"
#endif

namespace PHARE::solver
{
template<bool offload>
auto* constexpr solver_select()
{
    using CPU_SolverPPC_t = PHARE::solver::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;

#if defined(PHARE_WITH_GPU)
    using GPU_SolverPPC_t = PHARE::solver::gpu::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;

#elif defined(WITH_RAJA) and defined(WITH_UMPIRE)
    using GPU_SolverPPC_t = PHARE::solver::llnl::SolverPPC<HybridModel_t, PHARE::amr::SAMRAI_Types>;

#elif defined(WITH_RAJA) or defined(WITH_UMPIRE)
#error // invalid, both RAJA and UMPIRE are required together.
#endif

    return static_cast<std::conditional_t<offload, GPU_SolverPPC_t, CPU_SolverPPC_t>*>(nullptr);
}

template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_,
         bool offload = false>
struct PHARE_Types
{
    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    // core deps
    using core_types   = PHARE::core::PHARE_Types<dimension, interp_order, offload>;
    using VecField_t   = typename core_types::VecField_t;
    using Electromag_t = typename core_types::Electromag_t;
    using Ions_t       = typename core_types::Ions_t;
    using GridLayout_t = typename core_types::GridLayout_t;
    using Electrons_t  = typename core_types::Electrons_t;
    // core deps

    using IPhysicalModel = PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>;
    using HybridModel_t  = PHARE::solver::HybridModel<GridLayout_t, Electromag_t, Ions_t,
                                                     Electrons_t, PHARE::amr::SAMRAI_Types>;
    using MHDModel_t = PHARE::solver::MHDModel<GridLayout_t, VecField_t, PHARE::amr::SAMRAI_Types>;

    using SolverPPC_t = decltype(*solver_select<offload>());

    using SolverMHD_t = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types>;
    using LevelInitializerFactory_t = PHARE::solver::LevelInitializerFactory<HybridModel_t>;

    // amr deps
    using amr_types        = PHARE::amr::PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using RefinementParams = typename amr_types::RefinementParams;

    using MessengerFactory // = amr/solver bidirectional dependency
        = PHARE::amr::MessengerFactory<MHDModel_t, HybridModel_t, RefinementParams>;
    // amr deps

    using MultiPhysicsIntegrator
        = PHARE::solver::MultiPhysicsIntegrator<MessengerFactory, LevelInitializerFactory_t,
                                                PHARE::amr::SAMRAI_Types>;
};

} // namespace PHARE::solver

#endif // PHARE_SOLVER_INCLUDE_H
