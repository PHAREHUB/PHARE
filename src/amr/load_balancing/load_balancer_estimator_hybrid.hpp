#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP

#include "load_balancer_estimator.hpp"


namespace PHARE::amr
{
class LoadBalancerEstimatorHybrid : public LoadBalancerEstimator
{
public:
    // LoadBalancerEstimatorHybrid(){};

    ~LoadBalancerEstimatorHybrid() = default; // the implementation of a virtual class NEEDS a dtor

    std::string name() const override;
    //  void estimate() override;

    // private:
};


inline std::string LoadBalancerEstimatorHybrid::name() const
{
    return std::string("Hybrid");
}


// inline void LoadBalancerEstimatorHybrid::estimate() {}




} // namespace PHARE::amr

#endif
