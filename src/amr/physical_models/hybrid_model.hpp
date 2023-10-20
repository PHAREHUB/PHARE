#ifndef PHARE_HYBRID_MODEL_HPP
#define PHARE_HYBRID_MODEL_HPP

#include <string>

#include "initializer/data_provider.hpp"
#include "core/models/hybrid_state.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"
#include "amr/resources_manager/resources_manager.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "core/data/vecfield/vecfield.hpp"

namespace PHARE::solver
{
/**
 * @brief The HybridModel class is a concrete implementation of a IPhysicalModel. The class
 * holds a HybridState and a ResourcesManager.
 */
template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types>
class HybridModel : public IPhysicalModel<AMR_Types>
{
public:
    using type_list = PHARE::core::type_list<GridLayoutT, Electromag, Ions, Electrons, AMR_Types>;
    using Interface = IPhysicalModel<AMR_Types>;
    using amr_types = AMR_Types;
    using patch_t   = typename AMR_Types::patch_t;
    using level_t   = typename AMR_Types::level_t;
    static const std::string model_name;
    using gridlayout_type           = GridLayoutT;
    using electromag_type           = Electromag;
    using vecfield_type             = typename Electromag::vecfield_type;
    using field_type                = typename vecfield_type::field_type;
    using ions_type                 = Ions;
    using particle_array_type       = typename Ions::particle_array_type;
    using resources_manager_type    = amr::ResourcesManager<gridlayout_type>;
    static constexpr auto dimension = GridLayoutT::dimension;
    using ParticleInitializerFactory
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;


    core::HybridState<Electromag, Ions, Electrons> state;
    std::shared_ptr<resources_manager_type> resourcesManager;


    virtual void initialize(level_t& level) override;


    /**
     * @brief allocate uses the ResourcesManager to allocate HybridState physical quantities on
     * the given Patch at the given allocateTime
     */
    virtual void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
    }


    auto patch_data_ids() const { return resourcesManager->restart_patch_data_ids(*this); }


    /**
     * @brief fillMessengerInfo describes which variables of the model are to be initialized or
     * filled at ghost nodes.
     */
    virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    auto setOnPatch(patch_t& patch) { return resourcesManager->setOnPatch(patch, *this); }


    HybridModel(PHARE::initializer::PHAREDict const& dict,
                std::shared_ptr<resources_manager_type> const& _resourcesManager)
        : IPhysicalModel<AMR_Types>{model_name}
        , state{dict}
        , resourcesManager{std::move(_resourcesManager)}
    {
    }


    virtual ~HybridModel() override {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    bool isUsable() const { return state.isUsable(); }

    bool isSettable() const { return state.isSettable(); }

    auto getCompileTimeResourcesUserList() const { return std::forward_as_tuple(state); }

    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(state); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
};




//-------------------------------------------------------------------------
//                             definitions
//-------------------------------------------------------------------------


template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types>::initialize(level_t& level)
{
    for (auto& patch : level)
    {
        // first initialize the ions
        auto layout = amr::layoutFromPatch<gridlayout_type>(*patch);
        auto& ions  = state.ions;
        auto _      = this->resourcesManager->setOnPatch(*patch, state.electromag, state.ions);

        for (auto& pop : ions)
        {
            auto const& info         = pop.particleInitializerInfo();
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);
        }

        state.electromag.initialize(layout);
    }


    resourcesManager->registerForRestarts(*this);
}



template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    hybridInfo.modelMagnetic        = core::VecFieldNames{state.electromag.B};
    hybridInfo.modelElectric        = core::VecFieldNames{state.electromag.E};
    hybridInfo.modelIonDensity      = state.ions.densityName();
    hybridInfo.modelIonBulkVelocity = core::VecFieldNames{state.ions.velocity()};
    hybridInfo.modelCurrent         = core::VecFieldNames{state.J};

    hybridInfo.initElectric.emplace_back(core::VecFieldNames{state.electromag.E});
    hybridInfo.initMagnetic.emplace_back(core::VecFieldNames{state.electromag.B});

    hybridInfo.ghostElectric.push_back(hybridInfo.modelElectric);
    hybridInfo.ghostMagnetic.push_back(hybridInfo.modelMagnetic);
    hybridInfo.ghostCurrent.push_back(core::VecFieldNames{state.J});
    hybridInfo.ghostBulkVelocity.push_back(hybridInfo.modelIonBulkVelocity);


    auto transform_ = [](auto& ions, auto& inserter) {
        std::transform(std::begin(ions), std::end(ions), std::back_inserter(inserter),
                       [](auto const& pop) { return pop.name(); });
    };
    transform_(state.ions, hybridInfo.interiorParticles);
    transform_(state.ions, hybridInfo.levelGhostParticlesOld);
    transform_(state.ions, hybridInfo.levelGhostParticlesNew);
    transform_(state.ions, hybridInfo.patchGhostParticles);
}



template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types>
const std::string HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types>::model_name
    = "HybridModel";

template<typename... Args>
HybridModel<Args...> hybrid_model_from_type_list(core::type_list<Args...>);

template<typename TypeList>
struct type_list_to_hybrid_model
{
    using type = decltype(hybrid_model_from_type_list(std::declval<TypeList>()));
};

template<typename TypeList>
using type_list_to_hybrid_model_t = typename type_list_to_hybrid_model<TypeList>::type;

} // namespace PHARE::solver

#endif
