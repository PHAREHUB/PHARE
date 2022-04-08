
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
    static std::unique_ptr<IWorkLoadEstimator> create(std::string modelName)
    {
        if (modelName == "HybridModel")
        {
            auto wl = std::make_unique<HybridWorkLoadEstimator<PHARE_T>>();
            wl->set_strategy("NPPC");
            // TODO la methode "create" doit prendre un second arg qui est la strat_name
            return wl;
        }
        else
            return {};
    }
};
} // namespace PHARE::amr

#endif
