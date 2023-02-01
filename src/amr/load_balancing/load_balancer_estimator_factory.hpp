#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_FACTORY_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_FACTORY_HPP

#include <memory>
#include <string>

#include "load_balancer_estimator.hpp"
#include "load_balancer_estimator_hybrid.hpp"



namespace PHARE::amr
{
template<std::size_t dim>
class LoadBalancerEstimatorFactory
{
public:
    static std::unique_ptr<LoadBalancerEstimator<dim>> create(std::string modelName)
    {
        if (modelName == "HybridModel")
            return std::make_unique<LoadBalancerEstimatorHybrid<dim>>();

        return nullptr;
    };

private:
};

} // namespace PHARE::amr

#endif
