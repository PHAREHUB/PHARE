#ifndef PHARE_TYPES_HPP
#define PHARE_TYPES_HPP


#include "phare_solver.hpp"
#include "initializer/data_provider.hpp"
#include "diagnostic/diagnostics.hpp"
#include "restarts/restarts.hpp"

namespace PHARE
{
template<std::size_t dimension_, std::size_t interp_order_, std::size_t nbRefinedPart_>
struct PHARE_Types
{
    static auto constexpr dimension     = dimension_;
    static auto constexpr interp_order  = interp_order_;
    static auto constexpr nbRefinedPart = nbRefinedPart_;

    using core_types                      = PHARE::core::PHARE_Types<dimension, interp_order>;
    using Array_t                         = typename core_types::Array_t;
    using Grid_t                          = typename core_types::Grid_t;
    using Field_t                         = typename core_types::Field_t;
    using VecField_t                      = typename core_types::VecField_t;
    using Electromag_t                    = typename core_types::Electromag_t;
    using Ions_t                          = typename core_types::Ions_t;
    using YeeLayout_t                     = typename core_types::YeeLayout_t;
    using GridLayout_t                    = typename core_types::GridLayout_t;
    using ParticleArray_t                 = typename core_types::ParticleArray_t;
    using MaxwellianParticleInitializer_t = typename core_types::MaxwellianParticleInitializer_t;
    using IonPopulation_t                 = typename core_types::IonPopulation_t;
    using Electrons_t                     = typename core_types::Electrons_t;
    using ParticleInitializerFactory      = typename core_types::ParticleInitializerFactory;



    using amr_types        = PHARE::amr::PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using hierarchy_t      = typename amr_types::hierarchy_t;
    using Splitter         = typename amr_types::Splitter_t;
    using RefinementParams = typename amr_types::RefinementParams;




    using solver_types   = PHARE::solver::PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using IPhysicalModel = typename solver_types::IPhysicalModel;
    using HybridModel_t  = typename solver_types::HybridModel_t;
    using MHDModel_t     = typename solver_types::MHDModel_t;
    using SolverPPC_t    = typename solver_types::SolverPPC_t;
    /*using SolverMHD_t      = typename solver_types::SolverMHD_t;*/
    using MessengerFactory          = typename solver_types::MessengerFactory;
    using LevelInitializerFactory_t = typename solver_types::LevelInitializerFactory_t;
    using MultiPhysicsIntegrator    = typename solver_types::MultiPhysicsIntegrator;
};

} // namespace PHARE

#endif // PHARE_TYPES_HPP
