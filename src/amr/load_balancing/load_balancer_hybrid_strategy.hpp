
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
public:
    virtual void compute(SAMRAI::hier::PatchLevel& level,
                         PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>& model)
        = 0;
};


} // namespace PHARE::amr

#endif
