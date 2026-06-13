#ifndef LOAD_BALANCER_HYBRID_STRATEGY_FACTORY_HPP
#define LOAD_BALANCER_HYBRID_STRATEGY_FACTORY_HPP


#include <memory>

#include "amr/load_balancing/load_balancer_strategy.hpp"
#include "amr/load_balancing/concrete_load_balancer_strategy_nppc.hpp"
#include "amr/load_balancing/concrete_load_balancer_strategy_homogeneous.hpp"


namespace PHARE::amr
{
template<typename HybridModel>
class LoadBalancerHybridStrategyFactory
{
public:
    static std::unique_ptr<LoadBalancerStrategy<HybridModel>> create(std::string strat_name,
                                                                     int const id)
    {
        if (strat_name == "nppc")
        {
            return std::make_unique<ConcreteLoadBalancerStrategyNPPC<HybridModel>>(id);
        }

        else if (strat_name == "homogeneous")
        {
            return std::make_unique<ConcreteLoadBalancerStrategyHomogeneous<HybridModel>>(id);
        }

        return nullptr;
    }
};


} // namespace PHARE::amr


#endif
