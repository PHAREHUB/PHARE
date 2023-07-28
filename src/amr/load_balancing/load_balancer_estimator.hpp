#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HPP

#include <cstddef>
#include <string>

#include <SAMRAI/hier/PatchLevel.h>

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/physical_model.hpp"



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
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>& model)
        = 0;


protected:
    int const id_;
};

} // namespace PHARE::amr

#endif
