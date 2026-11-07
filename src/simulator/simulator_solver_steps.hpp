#ifndef PHARE_SIMULATOR_SIMULATOR_SOLVER_STEPS_HPP
#define PHARE_SIMULATOR_SIMULATOR_SOLVER_STEPS_HPP

// the hybrid solver field steps, split from simulator_solvers.hpp for compiler memory

#include "phare_solver.hpp"
#include "amr/solvers/solver_ppc_steps.hpp"

namespace PHARE
{

template<auto opts>
void emit_hybrid_solver_steps()
{
    if constexpr (has_hybrid_v<opts>)
    {
        using Solver_t = solver::PHARE_Types<opts>::Hybrid::Solver_t;
        solver::SolverPPCSteps<Solver_t>::emit();
    }
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_SOLVER_STEPS_HPP*/
