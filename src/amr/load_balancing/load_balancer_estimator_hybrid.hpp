#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP

#include "load_balancer_estimator.hpp"


namespace PHARE::amr
{
template<std::size_t dim>
class LoadBalancerEstimatorHybrid : public LoadBalancerEstimator<dim>
{
public:
    // LoadBalancerEstimatorHybrid(){};

    ~LoadBalancerEstimatorHybrid() = default; // the implementation of a virtual class NEEDS a dtor

    std::string modelName() const override;
    //  void estimate() override;

    // private:
};


template<std::size_t dim>
inline std::string LoadBalancerEstimatorHybrid<dim>::modelName() const
{
    return std::string("Hybrid");
}


// inline void LoadBalancerEstimatorHybrid::estimate() {}




} // namespace PHARE::amr

#endif
