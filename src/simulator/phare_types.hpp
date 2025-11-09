#ifndef PHARE_TYPES_HPP
#define PHARE_TYPES_HPP


#include "phare_solver.hpp"


namespace PHARE
{
template<SimOpts opts>
struct PHARE_Types
{
    static auto constexpr dimension     = opts.dimension;
    static auto constexpr interp_order  = opts.interp_order;
    static auto constexpr nbRefinedPart = opts.nbRefinedPart;

    using core_types                 = core::PHARE_Types<opts>;
    using Array_t                    = core_types::Array_t;
    using Grid_t                     = core_types::Grid_t;
    using Field_t                    = core_types::Field_t;
    using VecField_t                 = core_types::VecField_t;
    using Electromag_t               = core_types::Electromag_t;
    using Ions_t                     = core_types::Ions_t;
    using GridLayout_t               = core_types::GridLayout_t;
    using ParticleArray_t            = core_types::ParticleArray_t;
    using IonPopulation_t            = core_types::IonPopulation_t;
    using Electrons_t                = core_types::Electrons_t;
    using ParticleInitializerFactory = core_types::ParticleInitializerFactory_t;

    using amr_types        = amr::PHARE_Types<opts>;
    using hierarchy_t      = amr_types::hierarchy_t;
    using Splitter         = amr_types::Splitter_t;
    using RefinementParams = amr_types::RefinementParams;

    using solver_types              = solver::PHARE_Types<opts>;
    using IPhysicalModel            = solver_types::IPhysicalModel;
    using HybridModel_t             = solver_types::HybridModel_t;
    using MHDModel_t                = solver_types::MHDModel_t;
    using SolverPPC_t               = solver_types::SolverPPC_t;
    using SolverMHD_t               = solver_types::SolverMHD_t;
    using MessengerFactory          = solver_types::MessengerFactory;
    using LevelInitializerFactory_t = solver_types::LevelInitializerFactory_t;
    using MultiPhysicsIntegrator    = solver_types::MultiPhysicsIntegrator_t;
};

} // namespace PHARE

#endif // PHARE_TYPES_HPP
