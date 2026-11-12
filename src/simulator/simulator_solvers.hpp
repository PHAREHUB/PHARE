#ifndef PHARE_SIMULATOR_SIMULATOR_SOLVERS_HPP
#define PHARE_SIMULATOR_SIMULATOR_SOLVERS_HPP

#include "simulator/simulator_decl.hpp"

#include <memory>

namespace PHARE
{

template<auto opts>
void Simulator<opts>::hybrid_solver_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_hybrid_v<opts>)
    {
        using SolverPPC = PHARETypes::Hybrid::Solver_t;
        multiphysInteg_->registerAndInitSolver(
            maxMHDLevel_, maxLevelNumber_ - 1,
            std::make_unique<SolverPPC>(dict["simulation"]["algo"]));
    }
}


template<auto opts>
void Simulator<opts>::mhd_solver_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_mhd_v<opts>)
    {
        using SolverMHD = PHARETypes::MHD::Solver_t;
        multiphysInteg_->registerAndInitSolver(
            0, maxMHDLevel_ - 1, std::make_unique<SolverMHD>(dict["simulation"]["algo"]));
    }
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_SOLVERS_HPP*/
