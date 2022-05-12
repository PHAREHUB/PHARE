
#ifndef PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP
#define PHARE_CONCRETE_HYBRID_WORKLOAD_STRATEGY_NPPC_HPP


#include "hybrid_workload_strategy.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"



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
    virtual std::string name() override { return std::string("NPPC"); };
};



template<typename HybridModel>
void ConcreteHybridWorkLoadEstimatorStrategyNPPC<HybridModel>::estimate(
    double* workload_val, HybridModel const& hybrid_model,
    typename Super::gridlayout_type const& layout)
{
    // bool constexpr c_ordering = false;
    auto wl_view = core::NdArrayView<dimension, float*>(workload_val, layout.nbrCells());
    //= core::NdArrayView<dimension, int, int*, c_ordering>(workload_val, layout.nbrCells()); TODO :
    // this new version of NDArrayView once it will be merged... or need a rebase


    if constexpr (dimension == 1)
    {
        auto endCell = layout.nbrCells()[0] - 1;

        auto start = layout.physicalStartIndex(QtyCentering::dual, PHARE::core::Direction::X);
        auto end   = layout.physicalEndIndex(QtyCentering::dual, PHARE::core::Direction::X);

        for (auto iCell = 0u, ix = start; iCell <= endCell; ++ix, ++iCell)
        {
            (void)ix;

            // TODO : get the appropriate value with the cell_map
            // workload_val[iCell] = 1.0;

            // utiliser les NDArrayView pour transformer (i,j,k) en 1 indice 1D du double*
        }
    }
}

} // namespace PHARE::core

#endif
