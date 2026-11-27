
#ifndef PHARE_LOAD_BALANCER_HYBRID_STRATEGY_HPP
#define PHARE_LOAD_BALANCER_HYBRID_STRATEGY_HPP

#include "amr/physical_models/physical_model.hpp"



namespace PHARE::amr
{
template<typename Model>
class LoadBalancerStrategy
{
    using amr_types = Model::amr_types;
    using level_t   = amr_types::level_t;

public:
    virtual ~LoadBalancerStrategy() {}

    virtual void compute(level_t& level, PHARE::solver::IPhysicalModel<amr_types>& model) = 0;
};


} // namespace PHARE::amr

#endif
