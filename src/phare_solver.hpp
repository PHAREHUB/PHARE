
#ifndef PHARE_SOLVER_INCLUDE_HPP
#define PHARE_SOLVER_INCLUDE_HPP

#include "phare_amr.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver_ppc.hpp"
#include "amr/level_initializer/level_initializer.hpp"
#include "amr/level_initializer/level_initializer_factory.hpp"
#include "amr/multiphysics_integrator.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/physical_model.hpp"

namespace PHARE::solver
{
template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_,
         template<typename> typename MHDTimeStepper>
struct PHARE_Types
{
    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    // core deps
    using core_types = PHARE::core::PHARE_Types<dimension, interp_order>;

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
    using SolverMHD_t = PHARE::solver::SolverMHD<MHDModel_t, PHARE::amr::SAMRAI_Types,
                                                 MHDTimeStepper<MHDModel_t>>;
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

#endif // PHARE_SOLVER_INCLUDE_HPP
