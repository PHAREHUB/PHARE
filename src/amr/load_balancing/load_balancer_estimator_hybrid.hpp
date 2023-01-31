#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP

#include <memory>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellData.h>

#include "load_balancer_estimator.hpp"
// #include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/physical_model.hpp"
// #include "core/data/ndarray/ndarray_vector.hpp"
// #include "amr/resources_manager/amr_utils.hpp"
// #include "core/utilities/point/point.hpp"
#include "load_balancer_hybrid_strategy.hpp"
#include "load_balancer_hybrid_strategy_factory.hpp"




namespace PHARE::amr
{
template<typename PHARE_T>
class LoadBalancerEstimatorHybrid : public LoadBalancerEstimator
{
    using HybridModel     = typename PHARE_T::HybridModel_t;
    using gridlayout_type = typename HybridModel::gridlayout_type;



public:
    // LoadBalancerEstimatorHybrid(int const id)
    LoadBalancerEstimatorHybrid(std::string strategy_name, int const id)
        : LoadBalancerEstimator{id}
        , strat_{LoadBalancerHybridStrategyFactory<PHARE_T>::create(strategy_name, id)}
    {
    }

    ~LoadBalancerEstimatorHybrid() = default; // the implementation of a virtual class NEEDS a dtor

    void estimate(SAMRAI::hier::PatchLevel& level,
                  solver::IPhysicalModel<SAMRAI_Types>& model) override;


private:
    std::unique_ptr<LoadBalancerHybridStrategy<PHARE_T>> strat_;
};



template<typename PHARE_T>
inline void
LoadBalancerEstimatorHybrid<PHARE_T>::estimate(SAMRAI::hier::PatchLevel& level,
                                               solver::IPhysicalModel<SAMRAI_Types>& model)
{
    strat_->compute(level, model);

    // static auto constexpr dimension = HybridModel::dimension;
    // auto& hybridModel               = dynamic_cast<HybridModel&>(model);
    // auto& resourcesManager          = hybridModel.resourcesManager;
    // auto& ions                      = hybridModel.state.ions;

    // bool constexpr c_ordering = false;

    // for (auto& patch : level)
    // {
    //     auto const& layout = layoutFromPatch<gridlayout_type>(*patch);

    //     auto patch_data_lb
    //         =
    //         dynamic_cast<SAMRAI::pdat::CellData<double>*>(patch->getPatchData(this->id_).get());
    //     auto load_balancer_val = patch_data_lb->getPointer();

    //     const auto& box = patch->getBox();

    //     auto lb_view = core::NdArrayView<dimension, double, double*,
    //     c_ordering>(load_balancer_val,
    //                                                                              layout.nbrCells());

    //     auto _ = resourcesManager->setOnPatch(*patch, ions);

    //     // The view on "load_balancer_val" do not have any ghost cells, meaning that this patch
    //     // data lies on a box defind by the nbr of cells in the physical domain
    //     // As the nbr of ghost cells in the box associated to the load balancer patch data is
    //     // null, we hence loop on an AMR index (then needing to firstly get the AMR box)
    //     // to get the AMRLocatStart and AMRLocalEnd.
    //     // Then, we build "by hand" the local index for the "lb_view" considering that the nbr
    //     // of ghost cell for this patch data is null.
    //     auto amr_box = layout.AMRBox();

    //     // The lb_view is a CellData, meaning that it is dual as the index of an amr box
    //     auto amr_start = amr_box.lower;
    //     auto amr_end   = amr_box.upper;

    //     // nbr of cells in the physical domain
    //     auto nbrCells = layout.nbrCells();


    //     if constexpr (dimension == 1)
    //     {
    //         for (auto i_loc = 0, i_amr = amr_start[0]; i_loc < nbrCells[0]; ++i_loc, ++i_amr)
    //         {
    //             auto amr_point{core::Point{i_amr}};
    //             auto amr_index = amr_point.template toArray<int>();

    //             auto nbr = 0;

    //             for (auto& pop : ions)
    //             {
    //                 const auto& domainParticles = pop.domainParticles();

    //                 nbr += domainParticles.nbr_particles_in(amr_index);
    //             }

    //             lb_view(i_loc) = nbr;
    //         }
    //     }


    //     if constexpr (dimension == 2)
    //     {
    //         for (auto i_loc = 0, i_amr = amr_start[0]; i_loc < nbrCells[0]; ++i_loc, ++i_amr)
    //         {
    //             for (auto j_loc = 0, j_amr = amr_start[1]; j_loc < nbrCells[1]; ++j_loc, ++j_amr)
    //             {
    //                 auto amr_point{core::Point{i_amr, j_amr}};
    //                 auto amr_index = amr_point.template toArray<int>();

    //                 auto nbr = 0;

    //                 for (auto& pop : ions)
    //                 {
    //                     const auto& domainParticles = pop.domainParticles();

    //                     nbr += domainParticles.nbr_particles_in(amr_index);
    //                 }

    //                 lb_view(i_loc, j_loc) = nbr;
    //             }
    //         }
    //     }


    //     if constexpr (dimension == 3)
    //     {
    //         for (auto i_loc = 0, i_amr = amr_start[0]; i_loc < nbrCells[0]; ++i_loc, ++i_amr)
    //         {
    //             for (auto j_loc = 0, j_amr = amr_start[1]; j_loc < nbrCells[1]; ++j_loc, ++j_amr)
    //             {
    //                 for (auto k_loc = 0, k_amr = amr_start[2]; k_loc < nbrCells[2];
    //                      ++k_loc, ++k_amr)
    //                 {
    //                     auto amr_point{core::Point{i_amr, j_amr, k_amr}};
    //                     auto amr_index = amr_point.template toArray<int>();

    //                     auto nbr = 0;

    //                     for (auto& pop : ions)
    //                     {
    //                         const auto& domainParticles = pop.domainParticles();

    //                         nbr += domainParticles.nbr_particles_in(amr_index);
    //                     }

    //                     lb_view(i_loc, j_loc, k_loc) = nbr;
    //                 }
    //             }
    //         }
    //     }
    // }
}

} // namespace PHARE::amr

#endif
