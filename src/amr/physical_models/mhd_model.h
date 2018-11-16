


#ifndef PHARE_MHD_MODEL_H
#define PHARE_MHD_MODEL_H

#include <string>

#include "evolution/messengers/mhd_messenger_info.h"
#include "physical_models/physical_model.h"
#include "tools/resources_manager.h"

#include "hybrid/hybrid_quantities.h"

namespace PHARE
{
using MHDQuantity = HybridQuantity;


template<typename VecFieldT>
class MHDState : public PhysicalState
{
public:
    /*virtual void allocate(ResourcesManager const& manager, SAMRAI::hier::Patch& patch) override
    {
        manager.allocate(B, patch);
        manager.allocate(V, patch);
    }*/

    VecFieldT B{"B", MHDQuantity::Vector::B};
    VecFieldT V{"V", MHDQuantity::Vector::V};
};




template<typename GridLayoutT, typename VecFieldT>
class MHDModel : public IPhysicalModel
{
public:
    static const std::string model_name;
    using resources_manager_type = ResourcesManager<GridLayoutT>;


    explicit MHDModel(std::shared_ptr<resources_manager_type> resourcesManager)
        : IPhysicalModel{model_name}
        , resourcesManager{std::move(resourcesManager)}
    {
    }


    virtual std::unique_ptr<IMessengerInfo> messengerInfo() const override
    {
        return std::make_unique<MHDMessengerInfo>();
    }

    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state.B, patch, allocateTime);
        resourcesManager->allocate(state.V, patch, allocateTime);
    }



    virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const override {}

    virtual ~MHDModel() override = default;

    MHDState<VecFieldT> state;
    std::shared_ptr<resources_manager_type> resourcesManager;
};

template<typename GridLayoutT, typename VecFieldT>
const std::string MHDModel<GridLayoutT, VecFieldT>::model_name = "MHDModel";



} // namespace PHARE

#endif
