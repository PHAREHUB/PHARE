
#ifndef PHARE_LOAD_BALANCER_HYBRID_STRATEGY_HPP
#define PHARE_LOAD_BALANCER_HYBRID_STRATEGY_HPP

#include <SAMRAI/hier/PatchLevel.h>

#include "amr/physical_models/physical_model.hpp"
#include "amr/types/amr_types.hpp"



namespace PHARE::amr
{
template<typename PHARE_T>
class LoadBalancerHybridStrategy
{
    using HybridModel = typename PHARE_T::HybridModel_t;
    using amr_types   = typename HybridModel::amr_types;
    using level_t     = typename amr_types::level_t;

public:
    virtual ~LoadBalancerHybridStrategy() {}

    virtual void compute(level_t& level, PHARE::solver::IPhysicalModel<amr_types>& model) = 0;
};


} // namespace PHARE::amr

#endif
