#ifndef PHARE_HYBRID_MODEL_HPP
#define PHARE_HYBRID_MODEL_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"

#include "core/utilities/variants.hpp"
#include "initializer/data_provider.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"

#include <string>
#include <variant>

namespace PHARE::solver
{
/**
 * @brief The HybridModel class is a concrete implementation of a IPhysicalModel. The class
 * holds a HybridState and a ResourcesManager.
 */
template<typename Types>
class HybridModel : public IPhysicalModel<typename Types::amr_types>
{
public:
    using GridLayoutT        = Types::GridLayout_t;
    using Electromag         = Types::Electromag_t;
    using Ions               = Types::Ions_t;
    using Electrons          = Types::Electrons_t;
    using AMR_Types          = Types::amr_types;
    using Grid_t             = Types::Grid_t;
    using ResourcesManager_t = Types::ResourcesManager_t;

    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::options.interp_order;

    using Interface              = IPhysicalModel<AMR_Types>;
    using amr_types              = AMR_Types;
    using electrons_t            = Electrons;
    using patch_t                = AMR_Types::patch_t;
    using level_t                = AMR_Types::level_t;
    using physical_quantity_type = core::HybridQuantity;
    using gridlayout_type        = GridLayoutT;
    using electromag_type        = Electromag;
    using vecfield_type          = Electromag::vecfield_type;
    using tensorfield_type       = Ions::tensorfield_type;
    using field_type             = vecfield_type::field_type;
    using grid_type              = Grid_t;
    using ions_type              = Ions;
    using particle_array_type    = Ions::particle_array_type;
    using resources_manager_type = ResourcesManager_t;
    using ParticleInitializerFactory
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;

    using Resources = std::variant<field_type, vecfield_type, tensorfield_type>;

    static constexpr std::string_view model_type_name = "HybridModel";
    static inline std::string const model_name{model_type_name};


    void initialize(level_t& level) override;

    void registerResources() override { resourcesManager->registerResources(*this); }


    /**
     * @brief allocate uses the ResourcesManager to allocate HybridState physical quantities on
     * the given Patch at the given allocateTime
     */
    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(*this, patch, allocateTime);
    }


    /**
     * @brief fillMessengerInfo describes which variables of the model are to be initialized or
     * filled at ghost nodes.
     */
    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;


    NO_DISCARD auto setOnPatch(patch_t& patch)
    {
        return resourcesManager->setOnPatch(patch, *this);
    }


    HybridModel(PHARE::initializer::PHAREDict const& dict,
                std::shared_ptr<resources_manager_type> const& _resourcesManager)
        : IPhysicalModel<AMR_Types>{model_name}
        , state{dict}
        , resourcesManager{_resourcesManager}
    {
        registerTemporaryRequirements();
    }


    virtual ~HybridModel() override {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto& getRunTimeResourcesViewList() { return resources; }
    NO_DISCARD auto& getRunTimeResourcesViewList() const { return resources; }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(state); }
    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state); }

    void registerTemporaryRequirements(auto&&... quantities);
    auto& getTmp(auto&& as) { return core::get_from_variants(resources, as); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



    core::HybridState<Electromag, Ions, Electrons> state;
    std::shared_ptr<resources_manager_type> resourcesManager;
    std::vector<Resources> resources;
    std::array<std::size_t, std::variant_size_v<Resources>> n_per_type = {0};

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
};


template<typename Types>
void HybridModel<Types>::registerTemporaryRequirements(auto&&... quantities)
{
    // todo figure out at runtime
    std::array<std::size_t, std::variant_size_v<Resources>> _n_per_type{2, state.ions.size() + 1,
                                                                        state.ions.size() + 1};

    auto const namer = []<auto i>(auto n) -> std::string {
        return "tmp" + std::array<std::string, 3>{"Field", "VecField", "TensorField"}[i]
               + std::to_string(n);
    };

    core::for_N<std::variant_size_v<Resources>>([&](auto i) {
        while (n_per_type[i] < _n_per_type[i])
        {
            auto const name = namer.template operator()<i>(n_per_type[i]);
            if constexpr (i == 0)
                resources.emplace_back(field_type{name, core::HybridQuantity::all_primal_field});
            else if constexpr (i == 1)
                resources.emplace_back(vecfield_type{name, core::HybridQuantity::Vector::V});
            else if constexpr (i == 2)
                resources.emplace_back(tensorfield_type{name, core::HybridQuantity::Tensor::M});
            else
                static_assert(core::dependent_false_v<Types>);

            ++n_per_type[i];
        }
    });

    assert(resources.size() == 2 + ((state.ions.size() + 1) * 2));
}

//-------------------------------------------------------------------------
//                             definitions
//-------------------------------------------------------------------------


template<typename Types>
void HybridModel<Types>::initialize(level_t& level)
{
    for (auto& patch : level)
    {
        // first initialize the ions
        auto layout = amr::layoutFromPatch<gridlayout_type>(*patch);
        auto& ions  = state.ions;
        auto _ = this->resourcesManager->setOnPatch(*patch, state.electromag, state.ions, state.J);

        for (auto& pop : ions)
        {
            auto const& info         = pop.particleInitializerInfo();
            auto particleInitializer = ParticleInitializerFactory::create(info);
            particleInitializer->loadParticles(pop.domainParticles(), layout);
        }

        state.electromag.initialize(layout);
    }
}



template<typename Types>
void HybridModel<Types>::fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& hybridInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

    // only the charge density is registered to the messenger and not the ion mass
    // density. Reason is that mass density is only used to compute the
    // total bulk velocity which is already registered to the messenger
    hybridInfo.modelMagnetic        = state.electromag.B.name();
    hybridInfo.modelElectric        = state.electromag.E.name();
    hybridInfo.modelIonDensity      = state.ions.chargeDensityName();
    hybridInfo.modelIonBulkVelocity = state.ions.velocity().name();
    hybridInfo.modelCurrent         = state.J.name();

    hybridInfo.initElectric.emplace_back(state.electromag.E.name());
    hybridInfo.initMagnetic.emplace_back(state.electromag.B.name());

    hybridInfo.ghostElectric.push_back(hybridInfo.modelElectric);
    hybridInfo.ghostMagnetic.push_back(hybridInfo.modelMagnetic);
    hybridInfo.ghostCurrent.push_back(state.J.name());
    hybridInfo.ghostBulkVelocity.push_back(hybridInfo.modelIonBulkVelocity);

    auto transform_ = [](auto& ions, auto& inserter) {
        std::transform(std::begin(ions), std::end(ions), std::back_inserter(inserter),
                       [](auto const& pop) { return pop.name(); });
    };
    transform_(state.ions, hybridInfo.interiorParticles);
    transform_(state.ions, hybridInfo.levelGhostParticlesOld);
    transform_(state.ions, hybridInfo.levelGhostParticlesNew);
    transform_(state.ions, hybridInfo.patchGhostParticles);

    for (auto const& pop : state.ions)
    {
        hybridInfo.ghostFlux.emplace_back(pop.flux().name());
        hybridInfo.sumBorderFields.emplace_back(pop.particleDensity().name());
        hybridInfo.sumBorderFields.emplace_back(pop.chargeDensity().name());
    }

    hybridInfo.maxBorderFields.emplace_back(state.ions.massDensity().name());
    hybridInfo.maxBorderFields.emplace_back(state.ions.chargeDensity().name());
    hybridInfo.maxBorderVecFields.emplace_back(state.ions.velocity().name());
}




template<typename Model>
auto constexpr is_hybrid_model(Model* m) -> decltype(m->model_type_name, bool())
{
    return Model::model_type_name == "HybridModel";
}

template<typename... Args>
auto constexpr is_hybrid_model(Args...)
{
    return false;
}

template<typename Model>
auto constexpr is_hybrid_model_v = is_hybrid_model(static_cast<Model*>(nullptr));



} // namespace PHARE::solver

#endif
