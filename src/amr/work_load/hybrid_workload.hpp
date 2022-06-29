
#ifndef PHARE_HYBRID_WORKLOAD_HPP
#define PHARE_HYBRID_WORKLOAD_HPP


#include <SAMRAI/hier/PatchLevel.h>

#include "amr/work_load/workload.hpp"
#include "amr/work_load/workload_base.hpp"
#include "core/work_load/hybrid_workload_strategy_factory.hpp"
#include "core/work_load/hybrid_workload_strategy.hpp"
#include "amr/resources_manager/amr_utils.hpp"



namespace PHARE::amr
{
template<typename PHARE_T>
class HybridWorkLoadEstimator : public IWorkLoadEstimator, private WorkLoadEstimatorBase<PHARE_T>
{
    using HybridModel     = typename PHARE_T::HybridModel_t;
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    HybridWorkLoadEstimator()
        : WorkLoadEstimatorBase<PHARE_T>{"hybrid"}
    {
    }

    ~HybridWorkLoadEstimator()
    {
        assert(this->variableDatabase_->checkVariableExists(this->workLoadName_));
        this->variableDatabase_->removeVariable(this->workLoadName_);
    }

    virtual void estimate(SAMRAI::hier::PatchLevel, // double*,
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const&) override;

    virtual void set_strategy(std::string) override;

    virtual std::string name() const override;



private:
    std::unique_ptr<core::HybridWorkLoadEstimatorStrategy<HybridModel>> strat_;
};



template<typename PHARE_T>
void HybridWorkLoadEstimator<PHARE_T>::set_strategy(std::string stratName)
{
    strat_ = core::HybridWorkLoadStrategyFactory<PHARE_T>::create(stratName);
};



template<typename PHARE_T>
void HybridWorkLoadEstimator<PHARE_T>::estimate(
    SAMRAI::hier::PatchLevel level,
    PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const& model)
{
    auto& hybridModel = dynamic_cast<HybridModel const&>(model);

    for (auto& patch : level)
    {
        auto const& layout = PHARE::amr::layoutFromPatch<gridlayout_type>(*patch);
        auto pd
            = dynamic_cast<SAMRAI::pdat::CellData<double>*>(patch->getPatchData(this->id_).get());
        // auto workload_val = pd->getPointer(); // TODO this is where the problem is with tests 47,
        // 48, 48, 54, 55, 56, 57

        // strat_->estimate(workload_val, hybridModel, layout);
    }
};



template<typename PHARE_T>
std::string HybridWorkLoadEstimator<PHARE_T>::name() const
{
    return std::string("HybridWorkLoadEstimator_") + strat_->name();
}




} // namespace PHARE::amr

#endif
