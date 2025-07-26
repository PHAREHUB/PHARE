

#ifndef PHARE_CONCRERTE_LOAD_BALANCER_HYBRID_STRATEGY_HOMOGENEOUS_HPP
#define PHARE_CONCRERTE_LOAD_BALANCER_HYBRID_STRATEGY_HOMOGENEOUS_HPP


#include <SAMRAI/pdat/CellData.h>
#include <string>

#include "amr/load_balancing/load_balancer_strategy.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/types/amr_types.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"



namespace PHARE::amr
{
template<typename Model>
class ConcreteLoadBalancerStrategyHomogeneous : public LoadBalancerStrategy<Model>
{
public:
    using gridlayout_type = typename Model::gridlayout_type;
    using amr_types       = typename Model::amr_types;
    using level_t         = typename amr_types::level_t;
    using cell_data_t     = SAMRAI::pdat::CellData<double>;

    ConcreteLoadBalancerStrategyHomogeneous(int const id)
        : id_{id}
    {
    }

    void compute(level_t& level, PHARE::solver::IPhysicalModel<amr_types>& model) override;


private:
    int const id_;
};



template<typename Model>
void ConcreteLoadBalancerStrategyHomogeneous<Model>::compute(
    level_t& level, PHARE::solver::IPhysicalModel<amr_types>& model)
{
    for (auto& patch : level)
    {
        auto const& layout     = layoutFromPatch<gridlayout_type>(*patch);
        auto patch_data_lb     = dynamic_cast<cell_data_t*>(patch->getPatchData(this->id_).get());
        auto load_balancer_val = patch_data_lb->getPointer();

        // The view on "load_balancer_val" do not have any ghost cells, meaning that this patch
        // data lies on a box defind by the nbr of cells in the physical domain
        // As the nbr of ghost cells in the box associated to the load balancer patch data is
        // null, we hence loop on an AMR index (then needing to firstly get the AMR box)
        // to get the AMRLocatStart and AMRLocalEnd.
        // Then, we build "by hand" the local index for the "lb_view" considering that the nbr
        // of ghost cell for this patch data is null.

        // The lb_view is a CellData, meaning that it is dual as the index of an amr box

        std::fill(load_balancer_val, load_balancer_val + core::product(layout.nbrCells()), 1);
    }
}

} // namespace PHARE::amr

#endif
