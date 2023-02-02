#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HPP

#include <cstddef>
#include <string>


namespace PHARE::amr
{
template<std::size_t dim>
class LoadBalancerEstimator
{
public:
    // LoadBalancerEstimator() {}

    virtual ~LoadBalancerEstimator() = default;

    virtual std::string modelName() const = 0;
    // virtual void estimate()    = 0;

private:
};

} // namespace PHARE::amr

#endif
