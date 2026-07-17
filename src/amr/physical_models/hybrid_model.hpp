#ifndef PHARE_HYBRID_MODEL_HPP
#define PHARE_HYBRID_MODEL_HPP


#include "core/def.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/utilities/thread_pool.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"

#include "initializer/data_provider.hpp"

#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/resources_manager.hpp"


#include "hybrid_model_storage.hpp"

#include <string>
#include <stdexcept>
#include <type_traits>

namespace PHARE::solver
{

// for optional compile time per patch resources - eg grids for reducing during diags
template<typename AMR_Types, typename Ions, typename Grid_t>
class HybridModelBase
    : public IPhysicalModel<AMR_Types>,
      public hybrid_model_storage<Ions::particle_array_type::layout_mode, Ions, Grid_t>
{
public:
    using Super     = IPhysicalModel<AMR_Types>;
    using storage_t = hybrid_model_storage<Ions::particle_array_type::layout_mode, Ions, Grid_t>;
    using resources_manager_type = storage_t::resources_manager_type;

    HybridModelBase(std::string const& model_name)
        : Super{model_name}
        , storage_t{}
    {
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return storage_t::getCompileTimeResourcesViewList();
    }
    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return storage_t::getCompileTimeResourcesViewList();
    }
};


/**
 * @brief The HybridModel class is a concrete implementation of a IPhysicalModel. The class
 * holds a HybridState and a ResourcesManager.
 */
template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
class HybridModel : public HybridModelBase<AMR_Types, Ions, Grid_t>
{
public:
    static constexpr auto dimension = GridLayoutT::dimension;

    using Interface              = IPhysicalModel<AMR_Types>;
    using Super                  = HybridModelBase<AMR_Types, Ions, Grid_t>;
    using amr_types              = AMR_Types;
    using electrons_t            = Electrons;
    using patch_t                = AMR_Types::patch_t;
    using level_t                = AMR_Types::level_t;
    using State_t                = core::HybridState<Electromag, Ions, Electrons>;
    using physical_quantity_type = core::HybridQuantity;
    using gridlayout_type        = GridLayoutT;
    using electromag_type        = Electromag;
    using vecfield_type          = Electromag::vecfield_type;
    using field_type             = vecfield_type::field_type;
    using grid_type              = Grid_t;
    using ions_type              = Ions;
    using particle_array_type    = Ions::particle_array_type;
    using resources_manager_type = Super::resources_manager_type;
    using ParticleInitializerFactory_t
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;
    // *_rt = real type, not tiled!
    using Field_rt       = Super::storage_t::field_type;
    using VecField_rt    = Super::storage_t::vecfield_type;
    using TensorField_rt = Super::storage_t::tensorfield_type;

    static constexpr std::string_view model_type_name = "HybridModel";
    static inline std::string const model_name{model_type_name};



    State_t state;
    std::shared_ptr<resources_manager_type> resourcesManager;


    void initialize(level_t& level) override;


    /**
     * @brief allocate uses the ResourcesManager to allocate HybridState physical quantities on
     * the given Patch at the given allocateTime
     */
    virtual void allocate(patch_t& patch, double const allocateTime) override
    {
        core::apply(getCompileTimeResourcesViewList(),
                    [&](auto& el) { resourcesManager->allocate(el, patch, allocateTime); });
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
        : Super{model_name}
        , state{dict}
        , resourcesManager{_resourcesManager}
    {
    }


    virtual ~HybridModel() override {}

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        auto tup = std::tuple_cat(std::forward_as_tuple(state),
                                  Super::getCompileTimeResourcesViewList());
        core::for_N<std::tuple_size_v<decltype(tup)>>([&](auto i) {
            static_assert(std::is_reference_v<std::tuple_element_t<i, decltype(tup)>>);
        });
        return tup;
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        auto tup = std::tuple_cat(std::forward_as_tuple(state),
                                  Super::getCompileTimeResourcesViewList());
        core::for_N<std::tuple_size_v<decltype(tup)>>([&](auto i) {
            static_assert(std::is_reference_v<std::tuple_element_t<i, decltype(tup)>>);
        });
        return tup;
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


    Field_rt static inline tmpField{"phare_scratch_field", core::HybridQuantity::Scalar::Vx};
    VecField_rt static inline tmpVec{"phare_scratch_vec_field", core::HybridQuantity::Vector::V};
    TensorField_rt static inline tmpTensor{"phare_scratch_tensor_field",
                                           core::HybridQuantity::Tensor::M};
    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
};




//-------------------------------------------------------------------------
//                             definitions
//-------------------------------------------------------------------------


template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types, Grid_t>::initialize(
    level_t& level)
{
    auto static const init_mode = core::get_env_as("PHARE_MODEL_INIT_MODE", std::size_t{0});

    auto const init = [this](auto& patch, auto& electromag, auto& ions) {
        auto const& _      = this->resourcesManager->setOnPatch(*patch, electromag, ions);
        auto const& layout = amr::layoutFromPatch<gridlayout_type>(*patch);
        for (auto& pop : ions)
            ParticleInitializerFactory_t::create(pop.particleInitializerInfo())
                ->loadParticles(pop.domainParticles(), layout);

        electromag.initialize(layout);
    };

    // first initialize the ions

    // DEADLOCK with python functions IF NOT GIL-LESS!
    if (init_mode == 3) // round robin thread pools per patch, thread per tile
    {
        if constexpr (core::is_tiled(particle_array_type::layout_mode))
        {
            using TileParticleArray_t = particle_array_type::per_tile_particles;
            using TileParticleInitializerFactory_t
                = core::ParticleInitializerFactory<TileParticleArray_t, gridlayout_type>;

            for (auto& patch : resourcesManager->enumerate(level, state.ions))
            {
                auto const pool_idx = core::ThreadPool::INSTANCE().first_ready_idx();
                auto& pool          = core::ThreadPool::INSTANCE().get_pool(pool_idx);
                auto const& layout  = amr::layoutFromPatch<gridlayout_type>(*patch);
                for (auto& pop : state.ions)
                    for (auto [tyle, particles] : enumerate_tiles(pop.domainParticles()))
                        pool.detach_task([layout = layout, tile = &tyle, pop = &pop]() {
                            TileParticleInitializerFactory_t::create(pop->particleInitializerInfo())
                                ->loadParticles((*tile)(), layout.copy_as(**tile));
                        });
            }

            for (auto& patch : resourcesManager->enumerate(level, state.electromag))
                state.electromag.initialize(amr::layoutFromPatch<gridlayout_type>(*patch));

            core::ThreadPool::INSTANCE().wait();
        }
        else
        {
            throw std::runtime_error("Wrong layout");
        }
    }


    // DEADLOCK with python functions IF NOT GIL-LESS!
    else if (init_mode == 2) // round robin thread pools
    {
        auto& pool = core::ThreadPool::INSTANCE();
        for (auto& patch : level)
            pool.async([&, ions = state.ions, electromag = state.electromag,
                        patch = &*patch]() mutable { init(patch, electromag, ions); });
    }


    // DEADLOCK with python functions IF NOT GIL-LESS!
    else if (init_mode == 1) // single thread pool
    {
        auto& pool = core::ThreadPool::pool(0);
        for (auto& patch : level)
            pool.detach_task([&, ions = state.ions, electromag = state.electromag,
                              patch = patch]() mutable { init(patch, electromag, ions); });
        pool.wait();
    }


    else if (init_mode == 0) // default serial
    {
        for (auto& patch : level)
            init(patch, state.electromag, state.ions);
    }

    else
        throw std::runtime_error("no HybridModel::initialize impl");
}



template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types, Grid_t>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
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
