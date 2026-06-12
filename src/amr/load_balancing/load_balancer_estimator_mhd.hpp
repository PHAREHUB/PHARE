#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_MHD_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_MHD_HPP

#include <memory>
#include <SAMRAI/hier/PatchLevel.h>


#include "amr/load_balancing/concrete_load_balancer_strategy_homogeneous.hpp"
#include "amr/physical_models/physical_model.hpp"

#include "load_balancer_estimator.hpp"
#include "load_balancer_strategy.hpp"


namespace PHARE::amr
{
template<typename PHARE_T>
class LoadBalancerEstimatorMHD : public LoadBalancerEstimator
{
    using MHDModel  = typename PHARE_T::MHDModel_t;
    using amr_types = typename MHDModel::amr_types;
    using level_t   = typename amr_types::level_t;

public:
    LoadBalancerEstimatorMHD(int const id)
        : LoadBalancerEstimator{id}
        , strat_{std::make_unique<ConcreteLoadBalancerStrategyHomogeneous<MHDModel>>(id)}
    {
    }

    // the implementation of a virtual class NEEDS a dtor
    ~LoadBalancerEstimatorMHD() = default;

    void estimate(level_t& level, solver::IPhysicalModel<amr_types>& model) override
    {
        strat_->compute(level, model);
    }


private:
    std::unique_ptr<ConcreteLoadBalancerStrategyHomogeneous<MHDModel>> strat_;
};

} // namespace PHARE::amr

#endif
