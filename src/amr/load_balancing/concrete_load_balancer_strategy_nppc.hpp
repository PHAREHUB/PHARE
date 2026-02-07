

#ifndef PHARE_CONCRERTE_LOAD_BALANCER_HYBRID_STRATEGY_NPPC_HPP
#define PHARE_CONCRERTE_LOAD_BALANCER_HYBRID_STRATEGY_NPPC_HPP

#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellData.h>

#include <string>

#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "amr/types/amr_types.hpp"
#include "amr/load_balancing/load_balancer_strategy.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/amr_utils.hpp"




namespace PHARE::amr
{
template<typename Model>
class ConcreteLoadBalancerStrategyNPPC : public LoadBalancerStrategy<Model>
{
public:
    using gridlayout_type = typename Model::gridlayout_type;
    using amr_types       = typename Model::amr_types;
    using level_t         = typename amr_types::level_t;
    using cell_data_t     = SAMRAI::pdat::CellData<double>;

    ConcreteLoadBalancerStrategyNPPC(int const id)
        : id_{id}
    {
    }

    void compute(level_t& level, PHARE::solver::IPhysicalModel<amr_types>& model) override;


private:
    int const id_;
};



template<typename Model>
void ConcreteLoadBalancerStrategyNPPC<Model>::compute(
    level_t& level, PHARE::solver::IPhysicalModel<amr_types>& model)
{
    bool static constexpr c_ordering = false;
    auto static constexpr dimension  = Model::dimension;

    auto& concreteModel    = dynamic_cast<Model&>(model);
    auto& resourcesManager = concreteModel.resourcesManager;
    auto& ions             = concreteModel.state.ions;

    for (auto& patch : level)
    {
        auto const& layout     = layoutFromPatch<gridlayout_type>(*patch);
        auto patch_data_lb     = dynamic_cast<cell_data_t*>(patch->getPatchData(this->id_).get());
        auto load_balancer_val = patch_data_lb->getPointer();
        auto lb_view = core::make_array_view<c_ordering>(load_balancer_val, layout.nbrCells());
        auto _       = resourcesManager->setOnPatch(*patch, ions);

        // The view on "load_balancer_val" do not have any ghost cells, meaning that this patch
        // data lies on a box defind by the nbr of cells in the physical domain
        // As the nbr of ghost cells in the box associated to the load balancer patch data is
        // null, we hence loop on an AMR index (then needing to firstly get the AMR box)
        // to get the AMRLocatStart and AMRLocalEnd.
        // Then, we build "by hand" the local index for the "lb_view" considering that the nbr
        // of ghost cell for this patch data is null.

        // The lb_view is a CellData, meaning that it is dual as the index of an amr box
        core::Box<std::uint32_t, dimension> local_box{
            core::Point{core::ConstArray<std::uint32_t, dimension>()},
            core::Point{
                core::generate([](auto const& nCell) { return nCell - 1; }, layout.nbrCells())}};

        auto amr_iter = layout.AMRBox().begin();
        auto lcl_iter = local_box.begin();

        for (; lcl_iter != local_box.end(); ++amr_iter, ++lcl_iter)
            lb_view(*lcl_iter) = core::sum_from(ions, [&](auto const& pop) {
                return pop.domainParticles().nbr_particles_in((*amr_iter).toArray());
            });
    }

    // TODO here, we have the lb_view value correctly set on all patches. we also know the id_
    // so this is where we should call setWorkloadPatchDataIndex... which is a method of the
    // CascadePartitioner lb_view is a local container containing the data the LoadBalancerNanager
    // knows the id, as well as the LoadBalancerEstimator and the LoadBalancerEstimator is a
    // CascadePartitioner
}

} // namespace PHARE::amr

#endif
