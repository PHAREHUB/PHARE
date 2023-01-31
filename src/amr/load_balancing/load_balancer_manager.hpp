#ifndef PHARE_LOAD_BALANCER_MANAGER_HPP
#define PHARE_LOAD_BALANCER_MANAGER_HPP

#include <memory>
#include <vector>

#include "phare_core.hpp"
#include "load_balancer_estimator.hpp"



namespace PHARE::amr
{
class LoadBalancerManager
{
public:
    LoadBalancerManager() = default;
    void addLoadBalancerEstimator(std::unique_ptr<amr::LoadBalancerEstimator> lbe);
    int numOfEstimators() const;

private:
    std::vector<std::shared_ptr<amr::LoadBalancerEstimator>> loadBalancerEstimators_;
};



inline void
LoadBalancerManager::addLoadBalancerEstimator(std::unique_ptr<amr::LoadBalancerEstimator> lbe)
{
    if (core::notIn(lbe, loadBalancerEstimators_))
    {
        loadBalancerEstimators_.push_back(std::move(lbe));
    }
    else
    {
        throw std::runtime_error("load balancer estimator " + lbe->name() + " already registered");
    }
}



inline int LoadBalancerManager::numOfEstimators() const
{
    return loadBalancerEstimators_.size() - 1;
}

} // namespace PHARE::amr

#endif
