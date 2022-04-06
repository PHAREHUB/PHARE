
#ifndef PHARE_WORKLOAD_FACTORY_HPP
#define PHARE_WORKLOAD_FACTORY_HPP


#include "workload.hpp"
#include "hybrid_workload.hpp"


namespace PHARE::amr
{
class WorkLoadEstimatorFactory
{
    public :
        static std::shared_ptr<IWorkLoadEstimator> create(std::string modelName)
        {
            if (modelName == "HybridModel")
                return std::make_shared<HybridWorkLoadEstimator>();
            else
                return {};
        }
};
} // namespace PHARE::amr

#endif
