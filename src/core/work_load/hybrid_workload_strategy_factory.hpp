
#ifndef PHARE_HYBRID_WORK_LOAD_STRATEGY_FACTORY_HPP
#define PHARE_HYBRID_WORK_LOAD_STRATEGY_FACTORY_HPP


#include "hybrid_workload_strategy.hpp"
#include "concrete_hybrid_workload_strategy_NPPC.hpp"


class HybridWorkLoadStrategyFactory
{
    public :
        static std::unique_ptr<PHARE::core::HybridWorkLoadEstimatorStrategy> create(std::string stratName)
            {
                if (stratName == "NPPC")
                    return std::make_unique<ConcreteHybridWorkLoadEstimatorStrategyNPPC>();
                else
                    return {};
            };
};

#endif
