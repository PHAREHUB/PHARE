#ifndef PHARE_HYBRID_MODEL_HPP
#define PHARE_HYBRID_MODEL_HPP


#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/mpi_utils.hpp"
#include "core/models/hybrid_state.hpp"
#include "core/data/ions/particle_initializers/particle_initializer_factory.hpp"

#include "initializer/data_provider.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/resources_manager/resources_manager.hpp"

#include <string>
#include <sstream>
#include <algorithm>

namespace PHARE::solver
{
/**
 * @brief The HybridModel class is a concrete implementation of a IPhysicalModel. The class
 * holds a HybridState and a ResourcesManager.
 */
template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
class HybridModel : public IPhysicalModel<AMR_Types>
{
public:
    static constexpr auto dimension = GridLayoutT::dimension;

    using Interface              = IPhysicalModel<AMR_Types>;
    using amr_types              = AMR_Types;
    using electrons_t            = Electrons;
    using patch_t                = AMR_Types::patch_t;
    using level_t                = AMR_Types::level_t;
    using gridlayout_type        = GridLayoutT;
    using electromag_type        = Electromag;
    using vecfield_type          = Electromag::vecfield_type;
    using field_type             = vecfield_type::field_type;
    using grid_type              = Grid_t;
    using ions_type              = Ions;
    using particle_array_type    = Ions::particle_array_type;
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, grid_type>;
    using ParticleInitializerFactory
        = core::ParticleInitializerFactory<particle_array_type, gridlayout_type>;

    static constexpr std::string_view model_type_name = "HybridModel";
    static inline std::string const model_name{model_type_name};


    core::HybridState<Electromag, Ions, Electrons> state;
    std::shared_ptr<resources_manager_type> resourcesManager;


    void initialize(level_t& level) override;


    /**
     * @brief allocate uses the ResourcesManager to allocate HybridState physical quantities on
     * the given Patch at the given allocateTime
     */
    virtual void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
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
        , resourcesManager{std::move(_resourcesManager)}
    {
    }


    virtual ~HybridModel() override {}

    std::string summarize(auto& hierarchy);

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




//-------------------------------------------------------------------------
//                             definitions
//-------------------------------------------------------------------------

struct ElectromagMinMax
{
    std::size_t constexpr static COMPS = 3; // vecfield
    core::Point<double, COMPS> minE, maxE;
    core::Point<double, COMPS> minB, maxB;

    static auto GET(auto& rm, auto& em, auto& lvl)
    {
        ElectromagMinMax mm;
        auto& [minE, maxE, minB, maxB] = mm;
        for (auto const& _ : rm.enumerate(lvl, em))
            for (std::size_t i = 0; i < COMPS; ++i)
            {
                minE[i] = std::min(minE[i], *std::min_element(em.E[i].begin(), em.E[i].end()));
                maxE[i] = std::max(maxE[i], *std::max_element(em.E[i].begin(), em.E[i].end()));
                minB[i] = std::min(minB[i], *std::min_element(em.B[i].begin(), em.B[i].end()));
                maxB[i] = std::max(maxB[i], *std::max_element(em.B[i].begin(), em.B[i].end()));
            }

        return mm;
    }

    auto collect() const
    {
        ElectromagMinMax mm;
        for (std::size_t i = 0; i < COMPS; ++i)
        {
            mm.minE[i] = core::mpi::min_on_rank0(minE[i]);
            mm.maxE[i] = core::mpi::max_on_rank0(maxE[i]);
            mm.minB[i] = core::mpi::min_on_rank0(minB[i]);
            mm.maxB[i] = core::mpi::max_on_rank0(maxB[i]);
        }
        return mm;
    }
};


inline std::ostream& operator<<(std::ostream& out, ElectromagMinMax const& mm)
{
    out << "Emin(" << mm.minE.str() << "), Emax(" << mm.maxE.str() << "), Bmin(" << mm.minB.str()
        << "), Bmax(" << mm.maxB.str() << ")";

    return out;
}



template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
std::string
HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types, Grid_t>::summarize(auto& hierarchy)
{
    std::stringstream ss;

    std::size_t total = 0;

    ss << std::endl;
    amr::onLevels(hierarchy, [&](auto& lvl) mutable {
        auto const lcl
            = core::sum_from(this->resourcesManager->enumerate(lvl, state.ions), [&](auto const&) {
                  return core::sum_from(
                      state.ions, [](auto const& pop) { return pop.domainParticles().size(); });
              });

        auto const on_lvl = core::mpi::sum_on_rank_0(lcl);
        total += on_lvl;

        auto const minmax
            = ElectromagMinMax::GET(*this->resourcesManager, state.electromag, lvl).collect();

        if (core::mpi::rank() == 0)
            ss << "lvl:" << lvl.getLevelNumber() << " " << minmax << " parts(" << on_lvl << "), "
               << std::endl;
    });

    if (core::mpi::rank() == 0)
        ss << "tot:" << total;

    return ss.str();
}


template<typename GridLayoutT, typename Electromag, typename Ions, typename Electrons,
         typename AMR_Types, typename Grid_t>
void HybridModel<GridLayoutT, Electromag, Ions, Electrons, AMR_Types, Grid_t>::initialize(
    level_t& level)
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
