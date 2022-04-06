
#ifndef PHARE_HYBRID_WORKLOAD_HPP
#define PHARE_HYBRID_WORKLOAD_HPP

#include <SAMRAI/hier/PatchLevel.h>

#include "amr/types/amr_types.hpp"
#include "amr/work_load/workload.hpp"
#include "core/work_load/hybrid_workload_strategy_factory.hpp"


namespace PHARE::amr
{
template<typename HybridModel>
class HybridWorkLoadEstimator : public IWorkLoadEstimator
{
using gridlayout_type = typename HybridModel::gridlayout_type;

public :
    virtual void estimate(SAMRAI::hier::PatchLevel levels, double* wl, PHARE::solver::IPhysicalModel<amr_t> const& model)
    {
        auto& hybridModel = dynamic_cast<HybridModel&>(model);

        for (auto& p : levels)
        {
            auto layout = PHARE::amr::layoutFromPatch<gridlayout_type>(patch);

            auto pd = p.getPatchData(this->getID());
            (void)pd;



            // TODO



        }
    };

    virtual void set_strategy(std::string stratName)
    {
        strat_ = HybridWorkLoadStrategyFactory::create(stratName);
    };

private:
    std::unique_ptr<HybridWorkLoadEstimatorStrategy> strat_;

};

} // namespace PHARE::amr

#endif
