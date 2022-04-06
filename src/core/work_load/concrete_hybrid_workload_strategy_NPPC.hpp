
#ifndef PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP
#define PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP


#include "hybrid_workload_strategy.hpp"
#include "phare_solver.hpp"



namespace PHARE::core
{
class ConcreteHybridWorkLoadEstimatorStrategyNPPC : public HybridWorkLoadEstimatorStrategy
{

public :
    void estimate(double* wl, PHARE::solver::HybridModel_t const& hybrid_model) { } ;
};
} // namespace PHARE::core

#endif
