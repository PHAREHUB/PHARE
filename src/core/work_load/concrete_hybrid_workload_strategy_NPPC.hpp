
#ifndef PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP
#define PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP


#include "hybrid_workload_strategy.hpp"



namespace PHARE::core
{
template<typename HybridModel>
class ConcreteHybridWorkLoadEstimatorStrategyNPPC
    : public HybridWorkLoadEstimatorStrategy<HybridModel>
{
    static auto constexpr dimension = HybridModel::dimension;
    using Super                     = HybridWorkLoadEstimatorStrategy<HybridModel>;

public:
    virtual void estimate(double*, HybridModel const&, typename Super::gridlayout_type const&);
    virtual std::string name() override { return std::string("HybridWorkLoadEstimator_NPPC"); };
};



template<typename HybridModel>
void ConcreteHybridWorkLoadEstimatorStrategyNPPC<HybridModel>::estimate(
    double* workload_val, HybridModel const& hybrid_model,
    typename Super::gridlayout_type const& layout)
{
    if constexpr (dimension == 1)
    {
        auto endCell = layout.nbrCells()[0] - 1;

        auto start
            = layout.physicalStartIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);
        auto end
            = layout.physicalEndIndex(PHARE::core::QtyCentering::dual, PHARE::core::Direction::X);

        for (auto iCell = 0u, ix = start; iCell <= endCell; ++ix, ++iCell)
        {
            (void)ix;

            // TODO : get the appropriate value with the cell_map
            workload_val[iCell] = 1.0;
        }
    }
}

} // namespace PHARE::core

#endif
