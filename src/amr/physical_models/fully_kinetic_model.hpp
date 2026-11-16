#ifndef PHARE_FULLY_KINETIC_MODEL_HPP
#define PHARE_FULLY_KINETIC_MODEL_HPP


#include "core/def.hpp"
#include "core/models/fully_kinetic_state.hpp"
#include "core/boundary/boundary_manager.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"

#include "initializer/data_provider.hpp"

#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/resources_manager.hpp"


#include <string>
#include <string_view>

namespace PHARE::solver
{
/**
 * @brief FullyKineticModel is the IPhysicalModel of a plasma where all species, electrons
 * included, are particles. It exposes the same types and state members as HybridModel (minus the
 * fluid electrons) so that the hybrid messenger machinery can be reused on it.
 */
template<typename GridLayoutT, typename Electromag, typename Ions, typename AMR_Types,
         typename Grid_t>
class FullyKineticModel : public IPhysicalModel<AMR_Types>
{
public:
    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::options.interp_order;

    using Interface              = IPhysicalModel<AMR_Types>;
    using amr_types              = AMR_Types;
    using patch_t                = AMR_Types::patch_t;
    using level_t                = AMR_Types::level_t;
    using physical_quantity_type = core::HybridQuantity;
    using gridlayout_type        = GridLayoutT;
    using electromag_type        = Electromag;
    using vecfield_type          = Electromag::vecfield_type;
    using field_type             = vecfield_type::field_type;
    using grid_type              = Grid_t;
    using ions_type              = Ions;
    using particle_array_type    = Ions::particle_array_type;
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, grid_type>;
    using boundary_manager_type  = core::BoundaryManager<gridlayout_type, field_type>;
    using ParticleInitializerFactory
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;

    static constexpr std::string_view model_type_name = "FullyKineticModel";
    static inline std::string const model_name{model_type_name};


    core::FullyKineticState<Electromag, Ions> state;
    double const c_over_va; // speed of light over the Alfven speed
    std::shared_ptr<resources_manager_type> resourcesManager;
    std::shared_ptr<boundary_manager_type> boundaryManager;


    FullyKineticModel(PHARE::initializer::PHAREDict const& dict,
                      std::shared_ptr<resources_manager_type> const& _resourcesManager)
        : IPhysicalModel<AMR_Types>{model_name}
        , state{dict}
        , c_over_va{dict["c_over_va"].template to<double>()}
        , resourcesManager{_resourcesManager}
    {
        // periodic only: no particle boundary conditions yet, and the Esirkepov segment would
        // need to be folded at a physical boundary. The messenger still needs a manager.
        if (dict.contains("grid"))
            core::validatePhysicalBoundariesDeclared<dimension>(dict["grid"]);
        if (dict.contains("grid") && dict["grid"].contains("periodicities"))
            for (auto const direction : std::string_view{"xyz"}.substr(0, dimension))
                if (!dict["grid"]["periodicities"][std::string{direction}].template to<bool>())
                    throw std::runtime_error(
                        "FullyKineticModel: only periodic domains are supported");
        boundaryManager = std::make_shared<boundary_manager_type>(
            PHARE::initializer::PHAREDict{}, std::vector<core::HybridQuantity::Scalar>{},
            std::vector<core::HybridQuantity::Vector>{});
    }

    ~FullyKineticModel() override {}

    NO_DISCARD double c2() const { return c_over_va * c_over_va; }


    /**
     * @brief initialize loads the particles of all species and the magnetic field.
     * The electric field starts at zero, the level initializer then moves it half a step back.
     */
    void initialize(level_t& level) override
    {
        for (auto& patch : level)
        {
            auto layout = amr::layoutFromPatch<gridlayout_type>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, state.electromag, state.ions);

            for (auto& pop : state.ions)
            {
                auto const& info         = pop.particleInitializerInfo();
                auto particleInitializer = ParticleInitializerFactory::create(info, pop.name());
                particleInitializer->loadParticles(pop.domainParticles(), layout);
            }

            state.electromag.initialize(layout);
            state.electromag.E.zero();
        }
    }


    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
    }


    /**
     * @brief fillMessengerInfo declares the same quantities as the hybrid model, plus the
     * particle current that is summed on patch borders.
     */
    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override
    {
        auto& fkInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);
        fillParticleModelMessengerInfo(state, fkInfo);
        fkInfo.modelCurrent = state.J.name();
    }


    NO_DISCARD auto setOnPatch(patch_t& patch)
    {
        return resourcesManager->setOnPatch(patch, *this);
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(state); }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
};




template<typename Model>
auto constexpr is_fully_kinetic_model(Model* m) -> decltype(m->model_type_name, bool())
{
    return Model::model_type_name == "FullyKineticModel";
}

template<typename... Args>
auto constexpr is_fully_kinetic_model(Args...)
{
    return false;
}

template<typename Model>
auto constexpr is_fully_kinetic_model_v = is_fully_kinetic_model(static_cast<Model*>(nullptr));


// models whose species are particles, and that share the hybrid diagnostics and messengers
template<typename Model>
auto constexpr is_particle_model_v = is_hybrid_model_v<Model> or is_fully_kinetic_model_v<Model>;


} // namespace PHARE::solver

#endif
