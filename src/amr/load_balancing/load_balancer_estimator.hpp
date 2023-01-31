#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HPP

#include <string>


namespace PHARE::amr
{
class LoadBalancerEstimator
{
public:
    // LoadBalancerEstimator() {}

    virtual ~LoadBalancerEstimator() = default;

    virtual std::string name() const = 0;
    // virtual void estimate()    = 0;

private:
};

} // namespace PHARE::amr

#endif
