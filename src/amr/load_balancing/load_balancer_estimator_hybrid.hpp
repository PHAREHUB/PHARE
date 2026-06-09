#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP

#include <memory>
#include <SAMRAI/hier/PatchLevel.h>


#include "amr/physical_models/physical_model.hpp"

#include "load_balancer_estimator.hpp"
#include "load_balancer_strategy.hpp"
#include "load_balancer_hybrid_strategy_factory.hpp"


namespace PHARE::amr
{
template<typename PHARE_T>
class LoadBalancerEstimatorHybrid : public LoadBalancerEstimator
{
    using HybridModel = typename PHARE_T::HybridModel_t;
    using amr_types   = typename HybridModel::amr_types;
    using level_t     = typename amr_types::level_t;

public:
    LoadBalancerEstimatorHybrid(std::string strategy_name, int const id)
        : LoadBalancerEstimator{id}
        , strat_{LoadBalancerHybridStrategyFactory<HybridModel>::create(strategy_name, id)}
    {
    }

    // the implementation of a virtual class NEEDS a dtor
    ~LoadBalancerEstimatorHybrid() = default;

    void estimate(level_t& level, solver::IPhysicalModel<amr_types>& model) override
    {
        strat_->compute(level, model);
    }


private:
    std::unique_ptr<LoadBalancerStrategy<HybridModel>> strat_;
};

} // namespace PHARE::amr

#endif
