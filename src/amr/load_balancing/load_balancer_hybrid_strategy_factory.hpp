
#ifndef LOAD_BALANCER_HYBRID_STRATEGY_FACTORY_HPP
#define LOAD_BALANCER_HYBRID_STRATEGY_FACTORY_HPP


#include <memory>

#include "amr/load_balancing/load_balancer_hybrid_strategy.hpp"
#include "amr/load_balancing/concrete_load_balancer_hybrid_strategy_nppc.hpp"



namespace PHARE::amr
{
template<typename PHARE_T>
class LoadBalancerHybridStrategyFactory
{
public:
    static std::unique_ptr<LoadBalancerHybridStrategy<PHARE_T>> create(std::string strat_name,
                                                                       int const id)
    {
        if (strat_name == "nppc")
        {
            return std::make_unique<ConcreteLoadBalancerHybridStrategyNPPC<PHARE_T>>(id);
        }

        else if (strat_name == "homogeneous")
        {
            return std::make_unique<ConcreteLoadBalancerHybridStrategyNPPC<PHARE_T>>(id);
        }

        return nullptr;
    }
};

} // namespace PHARE::amr

#endif
