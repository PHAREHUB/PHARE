
#ifndef PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP
#define PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP


#include "hybrid_workload_strategy.hpp"



namespace PHARE::core
{
template<typename HybridModel>
class ConcreteHybridWorkLoadEstimatorStrategyNPPC : public HybridWorkLoadEstimatorStrategy<HybridModel>
{
public :
    virtual void estimate(double* wl, HybridModel const& hybrid_model) override
    {



        // TODO



    };

    virtual std::string name() override { return std::string("NPPC"); };
};
} // namespace PHARE::core

#endif
