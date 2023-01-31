#ifndef PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP
#define PHARE_LOAD_BALANCER_ESTIMATOR_HYBRID_HPP

#include <memory>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellData.h>

#include "load_balancer_estimator.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "amr/data/particles/particles_data.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"


namespace PHARE::amr
{
template<typename PHARE_T>
class LoadBalancerEstimatorHybrid : public LoadBalancerEstimator
{
    using HybridModel     = typename PHARE_T::HybridModel_t;
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    LoadBalancerEstimatorHybrid(int const id)
        : LoadBalancerEstimator{id}
    {
    }

    ~LoadBalancerEstimatorHybrid() = default; // the implementation of a virtual class NEEDS a dtor

    void estim(SAMRAI::hier::PatchLevel& level,
               PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const& model) override;
};


template<typename PHARE_T>
inline void LoadBalancerEstimatorHybrid<PHARE_T>::estim(
    SAMRAI::hier::PatchLevel& level,
    PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const& model)
{
    static auto constexpr dimension = HybridModel::dimension;
    auto& hybridModel               = dynamic_cast<HybridModel const&>(model);
    bool constexpr c_ordering       = false;

    for (auto& patch : level)
    {
        auto const& layout = PHARE::amr::layoutFromPatch<gridlayout_type>(*patch);
        auto patch_data_lb
            = dynamic_cast<SAMRAI::pdat::CellData<double>*>(patch->getPatchData(this->id_).get());
        auto load_balancer_val = patch_data_lb->getPointer();

        auto lb_view = core::NdArrayView<dimension, double, double*, c_ordering>(load_balancer_val,
                                                                                 layout.nbrCells());


        using partPack = typename amr::ParticlesData<core::ParticleArray<dimension>>;


        auto patch_data_part
            = std::dynamic_pointer_cast<amr::ParticlesData<core::ParticleArray<dimension>>>(
                patch->getPatchData(0 /* TODO how to get the proper id ? */));
        auto particle_array
            = patch_data_part->getPointer(); // TODO in this ptr, we should use "domainParticles +
                                             // patchGhostParticles" no ?
        // auto pa_view        = core::NdArrayView<dimension, partPack, partPack*, c_ordering>(
        // particle_array, layout.nbrCells()); // TODO uneeded no ?

        // TODO for a proper load balancing, we should include the ghost cell (as they involve load
        // of the mpi ranks)
        //
        auto nbrCells = layout.nbrCells();

        if constexpr (dimension == 2)
        {
            for (auto iLoad_x = 0u; iLoad_x < nbrCells[0]; ++iLoad_x)
            {
                for (auto iLoad_y = 0u; iLoad_y < nbrCells[1]; ++iLoad_y)
                {
                    lb_view(iLoad_x, iLoad_y) = 1.0; // TODO this value should come from the cellMap
                }
            }
        }
    }
}

} // namespace PHARE::amr

#endif
