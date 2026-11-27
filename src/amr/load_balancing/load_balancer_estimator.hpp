#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HPP

#include "amr/physical_models/physical_model.hpp"
#include "amr/types/amr_types.hpp"

#include <SAMRAI/hier/PatchLevel.h>



namespace PHARE::amr
{
class LoadBalancerEstimator
{
public:
    LoadBalancerEstimator(int const id)
        : id_{id}
    {
    }

    virtual ~LoadBalancerEstimator() = default;

    virtual void estimate(SAMRAI::hier::PatchLevel& level,
                          solver::IPhysicalModel<amr::SAMRAI_Types>& model)
        = 0;


protected:
    int const id_;
};

} // namespace PHARE::amr

#endif
