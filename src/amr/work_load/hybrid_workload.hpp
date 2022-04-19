
#ifndef PHARE_HYBRID_WORKLOAD_HPP
#define PHARE_HYBRID_WORKLOAD_HPP

#include <SAMRAI/hier/PatchLevel.h>

#include "amr/work_load/workload.hpp"
#include "core/work_load/hybrid_workload_strategy_factory.hpp"



namespace PHARE::amr
{
template<typename PHARE_T>
class HybridWorkLoadEstimator : public IWorkLoadEstimator<PHARE_T>
{
    using HybridModel     = typename PHARE_T::HybridModel_t;
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    virtual void estimate(SAMRAI::hier::PatchLevel, double*,
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const&) override;
    virtual void set_strategy(std::string) override;
    std::string name() const override;

private:
    std::unique_ptr<HybridWorkLoadEstimatorStrategy<HybridModel>> strat_;
};



template<typename PHARE_T>
void HybridWorkLoadEstimator<PHARE_T>::set_strategy(std::string stratName)
{
    strat_ = HybridWorkLoadStrategyFactory<PHARE_T>::create(stratName);
};



template<typename PHARE_T>
void HybridWorkLoadEstimator<PHARE_T>::estimate(
    SAMRAI::hier::PatchLevel levels, double* workload_val,
    PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const& model)
{
    auto& hybridModel = dynamic_cast<HybridModel const&>(model);

    for (auto& patch : levels)
    {
        auto const& layout = PHARE::amr::layoutFromPatch<gridlayout_type>(*patch);
        auto pd
            = dynamic_cast<SAMRAI::pdat::CellData<double>*>(patch->getPatchData(this->id_).get());
        auto workload_val = pd->getPointer();



        strat_->estimate(workload_val, hybridModel, layout);
    }
};



template<typename PHARE_T>
std::string HybridWorkLoadEstimator<PHARE_T>::name() const
{
    return std::string("HybridworkLoadEstimator_") + strat_->name();
}

} // namespace PHARE::amr

#endif
