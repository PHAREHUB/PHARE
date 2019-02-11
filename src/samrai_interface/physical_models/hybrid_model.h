


#ifndef PHARE_HYBRID_MODEL_H
#define PHARE_HYBRID_MODEL_H

#include <string>

#include "data_provider.h"
#include "evolution/messengers/hybrid_messenger_info.h"
#include "models/hybrid_state.h"
#include "physical_models/physical_model.h"
#include "tools/resources_manager.h"

namespace PHARE
{
namespace amr_interface
{
    /**
     * @brief The HybridModel class is a concrete implementation of a IPhysicalModel. The class
     * holds a HybridState and a ResourcesManager.
     */
    template<typename GridLayoutT, typename Electromag, typename Ions, typename IonsInitializer>
    class HybridModel : public IPhysicalModel
    {
    public:
        static const std::string model_name;
        using gridLayout_type           = GridLayoutT;
        using electromag_type           = Electromag;
        using vecfield_type             = typename Electromag::vecfield_type;
        using ions_type                 = Ions;
        using resources_manager_type    = ResourcesManager<gridLayout_type>;
        static constexpr auto dimension = GridLayoutT::dimension;


        //! Physical quantities associated with hybrid equations
        core::HybridState<Electromag, Ions, IonsInitializer> state;

        //! ResourcesManager used for interacting with SAMRAI databases, patchdata etc.
        std::shared_ptr<resources_manager_type> resourcesManager;


        /**
         * @brief This constructor uses the IonInitializer to build the Ions of the HybridState and
         * stores the ResourcesManager for allocating data on patches with allocate()
         */
        /*[[deprecated]] HybridModel(IonsInitializer ionsInitializer,
                                   std::shared_ptr<resources_manager_type> resourcesManager)
            : IPhysicalModel{model_name}
            , state{std::move(ionsInitializer)}
            , resourcesManager{std::move(resourcesManager)}
        {
        }*/


        HybridModel(PHARE::initializer::PHAREDict<dimension> dict,
                    std::shared_ptr<resources_manager_type> resourcesManager)
            : IPhysicalModel{model_name}
            , state{dict}
            , resourcesManager{std::move(resourcesManager)}
        {
        }


        /**
         * @brief allocate uses the ResourcesManager to allocate HybridState physical quantities on
         * the given Patch at the given allocateTime
         */
        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) override
        {
            // TODO MODEL could do resourcesManager->allocate(state);
            resourcesManager->allocate(state.electromag, patch, allocateTime);
            resourcesManager->allocate(state.ions, patch, allocateTime);
        }


        /**
         * @brief fillMessengerInfo describes which variables of the model are to be initialized or
         * filled at ghost nodes.
         */
        virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const override
        {
            auto& modelInfo = dynamic_cast<HybridMessengerInfo&>(*info);

            modelInfo.modelMagnetic   = VecFieldDescriptor{state.electromag.B};
            modelInfo.modelElectric   = VecFieldDescriptor{state.electromag.E};
            modelInfo.modelIonBulk    = VecFieldDescriptor{state.ions.velocity()};
            modelInfo.modelIonDensity = FieldDescriptor{state.ions.densityName()};

            modelInfo.initElectric.emplace_back(state.electromag.E);
            modelInfo.initMagnetic.emplace_back(state.electromag.B);
            modelInfo.initIonBulk.emplace_back(state.ions.velocity());
            modelInfo.initIonDensity.emplace_back(state.ions.densityName());

            modelInfo.ghostElectric.push_back(modelInfo.modelElectric);
            modelInfo.ghostMagnetic.push_back(modelInfo.modelMagnetic);

            std::transform(std::begin(state.ions), std::end(state.ions),
                           std::back_inserter(modelInfo.interiorParticles),
                           [](auto const& pop) { return pop.name(); });


            std::transform(std::begin(state.ions), std::end(state.ions),
                           std::back_inserter(modelInfo.coarseToFineParticles),
                           [](auto const& pop) { return pop.name(); });


            std::transform(std::begin(state.ions), std::end(state.ions),
                           std::back_inserter(modelInfo.ghostParticles),
                           [](auto const& pop) { return pop.name(); });
        }

        virtual ~HybridModel() override = default;

        //@TODO make it a resourcesUser
    };

    template<typename GridLayoutT, typename Electromag, typename Ions, typename IonsInitializer>
    const std::string HybridModel<GridLayoutT, Electromag, Ions, IonsInitializer>::model_name
        = "HybridModel";

} // namespace amr_interface

} // namespace PHARE

#endif
