
#ifndef PHARE_WORKLOAD_FACTORY_HPP
#define PHARE_WORKLOAD_FACTORY_HPP

#include <memory>

#include "workload.hpp"
#include "hybrid_workload.hpp"
#include "core/work_load/hybrid_workload_strategy_factory.hpp"


namespace PHARE::amr
{
template<typename PHARE_T>
class WorkLoadEstimatorFactory
{
public:
    static std::unique_ptr<IWorkLoadEstimator<PHARE_T::dimension>> create(std::string modelName)
    {
        return std::make_unique<HybridWorkLoadEstimator<PHARE_T>>();
        // std::cout << modelName << std::endl;
        // if (modelName == "HybridModel")
        // {
        //     return std::make_unique<HybridWorkLoadEstimator<PHARE_T>>();
        // }
        // else
        // {
        //     return nullptr;
        // }
    }
};
} // namespace PHARE::amr

#endif
