
#ifndef PHARE_HYBRID_WORKLOAD_HPP
#define PHARE_HYBRID_WORKLOAD_HPP

#include <SAMRAI/hier/PatchLevel.h>

#include "amr/work_load/workload.hpp"
#include "core/work_load/hybrid_workload_strategy_factory.hpp"



namespace PHARE::amr
{
template<typename PHARE_T>
class HybridWorkLoadEstimator : public IWorkLoadEstimator
{
    using HybridModel     = typename PHARE_T::HybridModel_t;
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    virtual void estimate(SAMRAI::hier::PatchLevel levels, double* wl,
                          PHARE::solver::IPhysicalModel<amr_t> const& model) override
    {
        auto& hybridModel = dynamic_cast<HybridModel const&>(model);

        for (auto& patch : levels)
        {
            auto layout = PHARE::amr::layoutFromPatch<gridlayout_type>(*patch);

            auto pd = patch->getPatchData(this->getID());
            (void)pd;



            // TODO
        }
    };

    virtual void set_strategy(std::string stratName) override
    {
        strat_ = HybridWorkLoadStrategyFactory<PHARE_T>::create(stratName);
    };

    std::string name() const override
    {
        if (strat_ == nullptr)
            std::cout << "ta mere en slip" << std::endl;
        return std::string("HybridworkLoadEstimator_") + strat_->name();
    }

private:
    std::unique_ptr<HybridWorkLoadEstimatorStrategy<HybridModel>> strat_;
};

} // namespace PHARE::amr

#endif
