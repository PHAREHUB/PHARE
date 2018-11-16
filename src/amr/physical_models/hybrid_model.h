


#ifndef PHARE_HYBRID_MODEL_H
#define PHARE_HYBRID_MODEL_H

#include <string>

#include "evolution/messengers/hybrid_messenger_info.h"
#include "physical_models/physical_model.h"
#include "tools/resources_manager.h"

namespace PHARE
{
template<typename Electromag, typename Ions, typename IonsInitializer>
class HybridState : public PhysicalState
{
public:
    HybridState(IonsInitializer ionsInitializer)
        : ions{std::move(ionsInitializer)}
    {
    }


    Electromag electromag{"EM"};
    Ions ions;

    static constexpr std::size_t dimension = Electromag::dimension;
};



template<typename GridLayoutT, typename Electromag, typename Ions, typename IonsInitializer>
class HybridModel : public IPhysicalModel
{
public:
    static const std::string model_name;
    using gridLayout_type        = GridLayoutT;
    using resources_manager_type = ResourcesManager<gridLayout_type>;


    HybridState<Electromag, Ions, IonsInitializer> state;
    std::shared_ptr<resources_manager_type> resourcesManager;

    HybridModel(IonsInitializer ionsInitializer,
                std::shared_ptr<resources_manager_type> resourcesManager)
        : IPhysicalModel{model_name}
        , state{std::move(ionsInitializer)}
        , resourcesManager{std::move(resourcesManager)}
    {
    }


    virtual std::unique_ptr<IMessengerInfo> messengerInfo() const override
    {
        auto info = std::make_unique<HybridMessengerInfo>();

        // info->electricGhostNames = state.electromag.E.getComponentNames();
        // info->magneticGhostNames = state.electromag.B.getComponentNames();


        // auto E_IDs = resourcesManager.getIDs(state.electromag.E);
        // auto B_IDs = resourcesManager.getIDs(state.electromag.B);


        for (auto& pop : state.ions)
        {
            // auto Pop_Rho_ID = resourcesManager.getIDs(pop, pop.densityName());
        }
        // info->electricGhost(E_ID_s);
        // info->magneticGhost(B_ID_s);

        // info->electricDomain(E_ID_s);
        // info->magneticDomain(B_ID_s);

        //

        return info;
    }

    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state.electromag.E, patch, allocateTime);
        resourcesManager->allocate(state.electromag.B, patch, allocateTime);
        resourcesManager->allocate(state.ions, patch, allocateTime);
    }


    virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const override
    {
        //
        auto& modelInfo = dynamic_cast<HybridMessengerInfo&>(*info);

        auto magneticComponentNames = extractNames(state.electromag.B);
        auto electricComponentNames = extractNames(state.electromag.E);

        modelInfo.modelMagnetic = {state.electromag.B.name(), magneticComponentNames[0],
                                   magneticComponentNames[1], magneticComponentNames[2]};

        modelInfo.modelElectric = {state.electromag.E.name(), electricComponentNames[0],
                                   electricComponentNames[1], electricComponentNames[2]};


        modelInfo.initElectric.push_back(modelInfo.modelElectric);
        modelInfo.initMagnetic.push_back(modelInfo.modelMagnetic);

        modelInfo.ghostElectric.push_back(modelInfo.modelElectric);
        modelInfo.ghostMagnetic.push_back(modelInfo.modelMagnetic);
    }

    virtual ~HybridModel() override = default;

    //@TODO make it a resourcesUser
};

template<typename GridLayoutT, typename Electromag, typename Ions, typename IonsInitializer>
const std::string HybridModel<GridLayoutT, Electromag, Ions, IonsInitializer>::model_name
    = "HybridModel";



} // namespace PHARE

#endif
