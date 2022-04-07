
#ifndef PHARE_HYBRID_WORKLOAD_STRATEGY_FACTORY_HPP
#define PHARE_HYBRID_WORKLOAD_STRATEGY_FACTORY_HPP


#include "hybrid_workload_strategy.hpp"
#include "concrete_hybrid_workload_strategy_NPPC.hpp"



namespace PHARE::core
{
template<typename PHARE_T>
class HybridWorkLoadStrategyFactory
{
using HybridModel = typename PHARE_T::HybridModel_t;

public :
    static std::unique_ptr<PHARE::core::HybridWorkLoadEstimatorStrategy<HybridModel>> create(std::string stratName)
        {
            using HybridModel = typename PHARE_T::HybridModel_t;
            if (stratName == "NPPC")
                return std::make_unique<ConcreteHybridWorkLoadEstimatorStrategyNPPC<HybridModel>>();
            else
                return {};
        };
};
} // namespace PHARE::core

#endif
