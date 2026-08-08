#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "mpi/mpi_utils.hpp"

#include "core/def.hpp"
#include "core/utilities/variants.hpp"
#include "core/models/quantities/mhd_quantities.hpp"

#include "amr/amr_constants.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/messengers/field_operate_transaction.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include "diagnostic/diagnostic_props.hpp"
#include "diagnostic/computers/diagnostic_computers.hpp"

#include <SAMRAI/xfer/RefineAlgorithm.h>

#include <functional>
#include <memory>
#include <unordered_map>
#include <utility>
#include <variant>

namespace PHARE::diagnostic
{
// Generic Template declaration, to override per Concrete model type
class IModelView
{
public:
    inline virtual ~IModelView();
};
IModelView::~IModelView() {}


template<typename GridLayout>
NO_DISCARD auto getPatchProperties(std::string const& patchID, GridLayout const& grid)
{
    PatchProperties dict;
    dict["origin"] = grid.origin().toVector();
    dict["nbrCells"]
        = core::Point<std::uint32_t, GridLayout::dimension>{grid.nbrCells()}.toVector();
    dict["lower"]    = grid.AMRBox().lower.toVector();
    dict["upper"]    = grid.AMRBox().upper.toVector();
    dict["mpi_rank"] = static_cast<std::size_t>(mpi::rank());
    return dict;
}

NO_DISCARD inline auto getEmptyPatchProperties(PatchProperties const& dict = {})
{
    dict["origin"]   = std::vector<double>{};
    dict["nbrCells"] = std::vector<std::uint32_t>{};
    dict["lower"]    = std::vector<int>{};
    dict["upper"]    = std::vector<int>{};
    dict["mpi_rank"] = std::size_t{0};
    return dict;
}


template<typename Derived, typename Hierarchy, typename Model>
class BaseModelView : public IModelView
{
public:
    using GridLayout        = Model::gridlayout_type;
    using VecField          = Model::vecfield_type;
    using ResMan            = Model::resources_manager_type;
    using Field             = Model::field_type;
    using UserTypes_t       = ResMan::template ResourceUserTypesFor<Model>;
    using VecFieldData_t    = UserTypes_t::template UserTensorField_t</*rank=*/1>::patch_data_type;
    using TensorFieldData_t = UserTypes_t::template UserTensorField_t</*rank=*/2>::patch_data_type;

    static constexpr auto dimension = Model::dimension;


    BaseModelView(Hierarchy& hierarchy, Model& model)
        : model_{model}
        , hierarchy_{hierarchy}
    {
    }


    template<typename Action>
    void onLevels(Action&& action, std::size_t const minlvl = 0,
                  std::size_t const maxlvl = amr::MAX_LEVEL_IDX)
    {
        amr::onLevels(hierarchy_, std::forward<Action>(action), minlvl, maxlvl);
    }


    template<typename OnLevel, typename OrMissing>
    void onLevels(OnLevel&& onLevel, OrMissing&& orMissing, std::size_t const minlvl,
                  std::size_t const maxlvl)
    {
        amr::onLevels(hierarchy_, std::forward<OnLevel>(onLevel),
                      std::forward<OrMissing>(orMissing), minlvl, maxlvl);
    }


    template<typename Action>
    void visitHierarchy(Action&& action, int minLevel = 0, int maxLevel = 0)
    {
        amr::visitHierarchy<GridLayout>(hierarchy_, *model_.resourcesManager,
                                        std::forward<Action>(action), minLevel, maxLevel, model_);
    }

    NO_DISCARD auto boundaryConditions() const { return hierarchy_.boundaryConditions(); }
    NO_DISCARD auto domainBox() const { return hierarchy_.domainBox(); }
    NO_DISCARD auto origin() const { return std::vector<double>(dimension, 0); }


    NO_DISCARD auto getPatchProperties(std::string patchID, GridLayout const& grid) const
    {
        return diagnostic::getPatchProperties(patchID, grid);
    }

    NO_DISCARD static auto getEmptyPatchProperties(PatchProperties dict = {})
    {
        return diagnostic::getEmptyPatchProperties(dict);
    }

    NO_DISCARD bool hasTagsVectorFor(int ilevel, std::string patch_id) const
    {
        auto key = std::to_string(ilevel) + "_" + patch_id;
        return model_.tags.count(key);
    }

    NO_DISCARD auto& getTagsVectorFor(int ilevel, std::string patch_id) const
    {
        auto key = std::to_string(ilevel) + "_" + patch_id;
        return model_.tags.at(key);
    }


protected:
    Model& model_;
    Hierarchy& hierarchy_;

private:
    Derived& derived() { return static_cast<Derived&>(*this); }
    Derived const& derived() const { return static_cast<Derived const&>(*this); }
};


/** \brief {weak PatchLevel identity, RefineSchedule} pair - shared entry type for the per-level
 * schedule caches below (MTAlgo::MTschedules, KineticEnergyFluxAlgo::KEFschedules).
 */
struct ScheduleEntry
{
    std::weak_ptr<SAMRAI::hier::PatchLevel> level; // invalidated if the Level is destroyed
    std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule;
};

/** \brief returns the RefineSchedule cached for level ilvl, rebuilding it if none exists yet or
 * if the PatchLevel at ilvl has been replaced since (regrid) - identity is checked via
 * ScheduleEntry::level rather than a bare .expired(), since the cached schedule itself likely
 * holds a strong ref to the old PatchLevel, which would keep expired() false forever.
 */
template<typename FieldData_t, typename PlusEqualsOp>
auto& getOrCreateBorderSumSchedule(auto& hierarchy, int const ilvl, auto& algo,
                                   std::map<int, ScheduleEntry>& schedules)
{
    auto const level   = hierarchy.getPatchLevel(ilvl);
    auto schedule_iter = schedules.find(ilvl);
    auto const create_schedule
        = schedule_iter == schedules.end() or schedule_iter->second.level.lock() != level;

    if (create_schedule)
        schedule_iter
            = schedules
                  .insert_or_assign(
                      ilvl, ScheduleEntry{level,
                                          algo->createSchedule(
                                              level, 0,
                                              std::make_shared<amr::FieldBorderOpTransactionFactory<
                                                  FieldData_t, PlusEqualsOp>>())})
                  .first;

    return *schedule_iter->second.schedule;
}


template<typename Hierarchy, typename Model>
class ModelView;


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
class ModelView<Hierarchy, Model>
    : public BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>
{
    using Super          = BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>;
    using Field          = Model::field_type;
    using VecField       = Model::vecfield_type;
    using VecFieldData_t = Super::VecFieldData_t;
    using TensorFieldT   = Model::ions_type::tensorfield_type;

public:
    using GridLayout             = Super::GridLayout;
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;

    struct FluidAccessors;
    struct FluidComputers;

    ModelView(Hierarchy& hierarchy, Model& model)
        : Super{hierarchy, model}
    {
        declareMomentumTensorAlgos();
        declareKineticEnergyFluxAlgos();
    }

    ModelView(ModelView const&)            = delete;
    ModelView(ModelView&&)                 = default;
    ModelView& operator=(ModelView const&) = delete;
    ModelView& operator=(ModelView&&)      = default;

    NO_DISCARD VecField& getB() const { return this->model_.state.electromag.B; }
    NO_DISCARD VecField& getE() const { return this->model_.state.electromag.E; }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }

    NO_DISCARD auto& model() const { return this->model_; }


    auto& tmpField(std::size_t idx = 0)
    {
        return this->model_.getTmp(
            Field{"tmpField" + std::to_string(idx), core::HybridQuantity::all_primal_field});
    }

    auto& tmpVecField(std::size_t idx = 0)
    {
        return this->model_.getTmp(
            VecField{"tmpVecField" + std::to_string(idx), core::HybridQuantity::Vector::V});
    }

    template<std::size_t rank = 2>
    auto& tmpTensorField(std::size_t idx = 0)
    {
        static_assert(rank > 0 and rank < 3);

        if constexpr (rank == 1)
            return tmpVecField(idx);
        else
            return this->model_.getTmp(TensorFieldT{"tmpTensorField" + std::to_string(idx),
                                                    core::HybridQuantity::Tensor::M});
    }

    void fillPopMomTensor(auto& lvl, auto const time, auto const popidx);


    NO_DISCARD std::vector<VecField*> getElectromagFields() const
    {
        return {&this->model_.state.electromag.B, &this->model_.state.electromag.E};
    }



    void fillPopKineticEnergyFluxVector(auto& lvl, auto const time, auto const popidx);


protected:
    struct MTAlgo;
    struct KineticEnergyFluxAlgo;

    void declareKineticEnergyFluxAlgos();
    void declareMomentumTensorAlgos();

    std::vector<MTAlgo> MTAlgos;
    std::vector<KineticEnergyFluxAlgo> kineticEnergyFluxAlgos;
    std::unique_ptr<FluidAccessors> fluidAccessors_;
    std::unique_ptr<FluidComputers> fluidComputers_;
};



template<typename Hierarchy, typename Model>
    requires solver::is_mhd_model_v<Model>
class ModelView<Hierarchy, Model>
    : public BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>
{
    using Field    = Model::field_type;
    using VecField = Model::vecfield_type;

public:
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;
    using BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>::BaseModelView;

    struct MHDAccessors;
    struct MHDComputers;

    NO_DISCARD auto& model() const { return this->model_; }

    NO_DISCARD const Field& getRho() const { return this->model_.state.rho; }

    NO_DISCARD const VecField& getRhoV() const { return this->model_.state.rhoV; }

    NO_DISCARD const VecField& getB() const { return this->model_.state.B; }

    NO_DISCARD const Field& getEtot() const { return this->model_.state.Etot; }

    NO_DISCARD const VecField& getE() const
    {
        throw std::runtime_error("E not currently available in MHD diagnostics");
    }

    // for setBuffer function in visitHierarchy
    NO_DISCARD Field& getRho() { return this->model_.state.rho; }

    NO_DISCARD VecField& getRhoV() { return this->model_.state.rhoV; }

    NO_DISCARD VecField& getB() { return this->model_.state.B; }

    NO_DISCARD Field& getEtot() { return this->model_.state.Etot; }

    NO_DISCARD VecField& getE()
    {
        throw std::runtime_error("E not currently available in MHD diagnostics");
    }

    auto& tmpField(std::size_t idx = 0)
    {
        return this->model_.getTmp(
            Field{"tmpField" + std::to_string(idx), core::MHDQuantity::all_primal_field});
    }

    auto& tmpVecField(std::size_t idx = 0)
    {
        return this->model_.getTmp(
            VecField{"tmpVecField" + std::to_string(idx), core::MHDQuantity::Vector::V});
    }

    template<std::size_t rank = 2>
    auto& tmpTensorField(std::size_t idx = 0)
    {
        static_assert(rank == 1);
        return tmpVecField(idx);
    }

private:
    std::unique_ptr<MHDAccessors> mhdAccessors_;
    std::unique_ptr<MHDComputers> mhdComputers_;
};



template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
struct ModelView<Hierarchy, Model>::MTAlgo
{
    auto& getOrCreateSchedule(auto& hierarchy, int const ilvl)
    {
        using PlusEqualsOp = core::PlusEquals<typename VecField::value_type>;
        return getOrCreateBorderSumSchedule<typename Super::TensorFieldData_t, PlusEqualsOp>(
            hierarchy, ilvl, MTalgo, MTschedules);
    }

    std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> MTalgo
        = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
    std::map<int, ScheduleEntry> MTschedules;
};


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
struct ModelView<Hierarchy, Model>::KineticEnergyFluxAlgo
{
    auto& getOrCreateSchedule(auto& hierarchy, int const ilvl)
    {
        using PlusEqualsOp = core::PlusEquals<typename VecField::value_type>;
        return getOrCreateBorderSumSchedule<VecFieldData_t, PlusEqualsOp>(hierarchy, ilvl, KEFalgo,
                                                                          KEFschedules);
    }

    std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> KEFalgo
        = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
    std::map<int, ScheduleEntry> KEFschedules;
};


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
void ModelView<Hierarchy, Model>::declareMomentumTensorAlgos()
{
    auto& rm = *this->model_.resourcesManager;

    auto const dst_name = tmpTensorField(this->model_.state.ions.size()).name();

    for (std::size_t i = 0; i < this->model_.state.ions.size(); ++i)
    {
        auto& MTAlgo        = MTAlgos.emplace_back();
        auto const src_name = tmpTensorField(i).name();

        auto&& [idDst, idSrc] = rm.getIDsList(dst_name, src_name);
        MTAlgo.MTalgo->registerRefine(
            idDst, idSrc, idDst, nullptr,
            std::make_shared<amr::TensorFieldGhostInterpOverlapFillPattern<
                typename Super::GridLayout, /*rank_=*/2>>());
    }

    // can't create schedules here as the hierarchy has no levels yet
}


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
void ModelView<Hierarchy, Model>::declareKineticEnergyFluxAlgos()
{
    auto& rm            = *this->model_.resourcesManager;
    auto const dst_name = tmpVecField(this->model_.state.ions.size()).name();
    for (std::size_t i = 0; i < this->model_.state.ions.size(); ++i)
    {
        auto& kineticEnergyFluxAlgo = kineticEnergyFluxAlgos.emplace_back();
        auto const src_name         = tmpVecField(i).name();

        auto&& [idDst, idSrc] = rm.getIDsList(dst_name, src_name);
        kineticEnergyFluxAlgo.KEFalgo->registerRefine(
            idDst, idSrc, idDst, nullptr,
            std::make_shared<
                amr::TensorFieldGhostInterpOverlapFillPattern<GridLayout, /*rank_=*/1>>());
    }
    // can't create schedules here as the hierarchy has no levels yet
}


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
void ModelView<Hierarchy, Model>::fillPopMomTensor(auto& lvl, auto const time, auto const popidx)
{
    using value_type = TensorFieldT::value_type;
    auto constexpr N = core::detail::tensor_field_dim_from_rank<2>();

    auto& rm      = *this->model_.resourcesManager;
    auto& ions    = this->model_.state.ions;
    auto& src     = tmpTensorField(popidx);
    auto& scratch = tmpTensorField(ions.size());

    for (auto patch : rm.enumerate(lvl, src, scratch))
        for (std::uint8_t c = 0; c < N; ++c)
            std::memcpy(scratch[c].data(), src[c].data(), src[c].size() * sizeof(value_type));

    MTAlgos[popidx].getOrCreateSchedule(this->hierarchy_, lvl.getLevelNumber()).fillData(time);

    for (auto patch : rm.enumerate(lvl, src, scratch))
        for (std::uint8_t c = 0; c < N; ++c)
            std::memcpy(src[c].data(), scratch[c].data(), src[c].size() * sizeof(value_type));
}


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
void ModelView<Hierarchy, Model>::fillPopKineticEnergyFluxVector(auto& lvl, auto const time,
                                                                 auto const popidx)
{
    using value_type = TensorFieldT::value_type;
    auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();

    auto& rm      = *this->model_.resourcesManager;
    auto& ions    = this->model_.state.ions;
    auto& src     = tmpVecField(popidx);
    auto& scratch = tmpVecField(ions.size());

    for (auto patch : rm.enumerate(lvl, src, scratch))
        for (std::uint8_t c = 0; c < N; ++c)
            std::memcpy(scratch[c].data(), src[c].data(), src[c].size() * sizeof(value_type));

    kineticEnergyFluxAlgos[popidx]
        .getOrCreateSchedule(this->hierarchy_, lvl.getLevelNumber())
        .fillData(time);

    for (auto patch : rm.enumerate(lvl, src, scratch))
        for (std::uint8_t c = 0; c < N; ++c)
            std::memcpy(src[c].data(), scratch[c].data(), src[c].size() * sizeof(value_type));
}


/** \brief {name, accessor} registry for hybrid fluid quantities ("/ions/density",
 * "/ions/pop/<name>/momentum_tensor", ...), keyed by the full diagnostic.quantity path. Built
 * once per ModelView from the real population names/indices, so callers do a single map lookup
 * instead of looping populations to find which one's tree prefix matches - shared across writer
 * types (h5, vtk) since it only depends on Model_t, not on any writer.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
struct ModelView<Hierarchy, Model>::FluidAccessors
{
    using Field_t       = Model_t::field_type;
    using VecField_t    = Model_t::vecfield_type;
    using TensorField_t = Model_t::tensorfield_type;
    using FieldVariant  = std::variant<Field_t*, VecField_t*, TensorField_t*>;
    using FieldFn       = std::function<FieldVariant(Model_t&)>;

    struct Entry
    {
        std::string ownerKey; // patchAttributes bucket: "ion" or "fluid_<pop name>"
        std::string name;     // bare quantity name, e.g. "density", "momentum_tensor"
        FieldFn get;
    };

    static FluidAccessors& getOrCreateFor(auto& modelView)
    {
        if (!modelView.fluidAccessors_)
            modelView.fluidAccessors_ = std::make_unique<FluidAccessors>(modelView);
        return *modelView.fluidAccessors_;
    }

    FluidAccessors(auto& modelView);

    // action(field, name, ownerKey) - false if qty has no registered accessor
    template<typename Action>
    bool dispatch(std::string const& qty, Model_t& model, Action&& action) const
    {
        if (auto it = fields.find(qty); it != fields.end())
        {
            auto const& entry = it->second;
            std::visit([&](auto* field) { action(*field, entry.name, entry.ownerKey); },
                       entry.get(model));
            return true;
        }
        return false;
    }

    std::unordered_map<std::string, Entry> fields;
};


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
ModelView<Hierarchy, Model>::FluidAccessors::FluidAccessors(auto& modelView)
{
    auto& ions = modelView.getIons();

    fields["/ions/charge_density"] = {"ion", "charge_density", [](Model_t& model) -> FieldVariant {
                                          return &model.state.ions.chargeDensity();
                                      }};
    fields["/ions/mass_density"]   = {"ion", "mass_density", [](Model_t& model) -> FieldVariant {
                                        return &model.state.ions.massDensity();
                                    }};
    fields["/ions/bulkVelocity"]   = {"ion", "bulkVelocity", [](Model_t& model) -> FieldVariant {
                                        return &model.state.ions.velocity();
                                    }};
    fields["/ions/momentum_tensor"]
        = {"ion", "momentum_tensor", [n = ions.size()](Model_t& model) -> FieldVariant {
               TensorField_t key{"tmpTensorField" + std::to_string(n),
                                 core::HybridQuantity::Tensor::M};
               return &core::get_from_variants(model.getRunTimeResourcesViewList(), key);
           }};
    fields["/ions/kinetic_energy_flux_vector"]
        = {"ion", "kinetic_energy_flux_vector", [n = ions.size()](Model_t& model) -> FieldVariant {
               VecField_t key{"tmpVecField" + std::to_string(n), core::HybridQuantity::Vector::V};
               return &core::get_from_variants(model.getRunTimeResourcesViewList(), key);
           }};

    for (std::size_t idx = 0; idx < ions.size(); ++idx)
    {
        std::string const tree{"/ions/pop/" + ions[idx].name() + "/"};
        std::string const ownerKey{"fluid_" + ions[idx].name()};

        fields[tree + "density"] = {ownerKey, "density", [idx](Model_t& model) -> FieldVariant {
                                        return &model.state.ions[idx].particleDensity();
                                    }};
        fields[tree + "charge_density"]
            = {ownerKey, "charge_density", [idx](Model_t& model) -> FieldVariant {
                   return &model.state.ions[idx].chargeDensity();
               }};
        fields[tree + "flux"] = {ownerKey, "flux", [idx](Model_t& model) -> FieldVariant {
                                     return &model.state.ions[idx].flux();
                                 }};
        fields[tree + "momentum_tensor"]
            = {ownerKey, "momentum_tensor", [idx](Model_t& model) -> FieldVariant {
                   TensorField_t key{"tmpTensorField" + std::to_string(idx),
                                     core::HybridQuantity::Tensor::M};
                   return &core::get_from_variants(model.getRunTimeResourcesViewList(), key);
               }};
        fields[tree + "kinetic_energy_flux_vector"]
            = {ownerKey, "kinetic_energy_flux_vector", [idx](Model_t& model) -> FieldVariant {
                   VecField_t key{"tmpVecField" + std::to_string(idx),
                                  core::HybridQuantity::Vector::V};
                   return &core::get_from_variants(model.getRunTimeResourcesViewList(), key);
               }};
    }
}


/** \brief {name, computer} registry for hybrid fluid quantities that require an actual
 * computation pass ("/ions/momentum_tensor", "/ions/pop/<name>/kinetic_energy_flux_vector",
 * ...), keyed by the full diagnostic.quantity path. Built once per ModelView from the real
 * population names/references, so callers do a single map lookup instead of looping
 * populations to find which one's tree prefix matches. Writer-agnostic (takes minLvl/maxLvl/
 * timestamp as plain arguments rather than a concrete writer type) so it's shared across
 * writer types the same way FluidAccessors is.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
struct ModelView<Hierarchy, Model>::FluidComputers
{
    using Computer = std::function<void(ModelView&, std::size_t, std::size_t, double)>;

    static FluidComputers& getOrCreateFor(auto& modelView)
    {
        if (!modelView.fluidComputers_)
            modelView.fluidComputers_ = std::make_unique<FluidComputers>(modelView);
        return *modelView.fluidComputers_;
    }

    FluidComputers(auto& modelView);

    // returns false if qty has no registered computer
    bool compute(std::string const& qty, ModelView& modelView, std::size_t minLvl,
                 std::size_t maxLvl, double timestamp) const
    {
        auto it = computers.find(qty);
        if (it == computers.end())
            return false;
        it->second(modelView, minLvl, maxLvl, timestamp);
        return true;
    }

    std::unordered_map<std::string, Computer> computers;
};


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
ModelView<Hierarchy, Model>::FluidComputers::FluidComputers(auto& modelView)
{
    computers["/ions/momentum_tensor"]
        = [](ModelView& mv, std::size_t minLvl, std::size_t maxLvl, double timestamp) {
              compute_momentum_tensor(mv, minLvl, maxLvl, timestamp);
          };
    computers["/ions/kinetic_energy_flux_vector"]
        = [](ModelView& mv, std::size_t minLvl, std::size_t maxLvl, double timestamp) {
              compute_kinetic_energy_flux_vector(mv, minLvl, maxLvl, timestamp);
          };

    for (auto& pop : modelView.getIons())
    {
        computers["/ions/pop/" + pop.name() + "/momentum_tensor"]
            = [&pop](ModelView& mv, std::size_t minLvl, std::size_t maxLvl, double timestamp) {
                  compute_pop_momentum_tensor(mv, pop, minLvl, maxLvl, timestamp);
              };
        computers["/ions/pop/" + pop.name() + "/kinetic_energy_flux_vector"]
            = [&pop](ModelView& mv, std::size_t minLvl, std::size_t maxLvl, double timestamp) {
                  compute_pop_kinetic_energy_flux_vector(mv, pop, minLvl, maxLvl, timestamp);
              };
    }
}


/** \brief {name, accessor} registry for MHD quantities ("/mhd/rho", "/mhd/V", ...), keyed by the
 * full diagnostic.quantity path - mirrors ModelView<Hierarchy, HybridModel>::FluidAccessors, but
 * the MHD state has no per-population tree, so there is a single "mhd" bucket and one entry per
 * state field.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_mhd_model_v<Model>
struct ModelView<Hierarchy, Model>::MHDAccessors
{
    using Field_t      = Model_t::field_type;
    using VecField_t   = Model_t::vecfield_type;
    using FieldVariant = std::variant<Field_t*, VecField_t*>;
    using FieldFn      = std::function<FieldVariant(Model_t&)>;

    struct Entry
    {
        std::string ownerKey; // patchAttributes bucket - always "mhd"
        std::string name;     // bare quantity name, e.g. "rho", "V"
        FieldFn get;
    };

    static MHDAccessors& getOrCreateFor(auto& modelView)
    {
        if (!modelView.mhdAccessors_)
            modelView.mhdAccessors_ = std::make_unique<MHDAccessors>(modelView);
        return *modelView.mhdAccessors_;
    }

    // defined inline (rather than declared here + defined out-of-line further down, like
    // FluidAccessors does) - by the time this second ModelView partial specialization is parsed,
    // both the hybrid and mhd specializations exist in the TU, and clang's qualified-lookup for
    // an out-of-line nested-ctor definition (ModelView<H,M>::MHDAccessors::MHDAccessors(...))
    // fails to disambiguate an unqualified `ModelView` parameter type in that position (it did
    // work for FluidAccessors, parsed earlier, while only the hybrid specialization existed).
    // Defining it here avoids the qualified-name path entirely.
    MHDAccessors(auto&)
    {
        std::string const tree{"/mhd/"};

        fields[tree + "rho"]
            = {"mhd", "rho", [](Model_t& model) -> FieldVariant { return &model.state.rho; }};
        fields[tree + "V"]
            = {"mhd", "V", [](Model_t& model) -> FieldVariant { return &model.state.V; }};
        fields[tree + "P"]
            = {"mhd", "P", [](Model_t& model) -> FieldVariant { return &model.state.P; }};
        fields[tree + "rhoV"]
            = {"mhd", "rhoV", [](Model_t& model) -> FieldVariant { return &model.state.rhoV; }};
        fields[tree + "Etot"]
            = {"mhd", "Etot", [](Model_t& model) -> FieldVariant { return &model.state.Etot; }};
    }

    // action(field, name, ownerKey) - false if qty has no registered accessor
    template<typename Action>
    bool dispatch(std::string const& qty, Model_t& model, Action&& action) const
    {
        if (auto it = fields.find(qty); it != fields.end())
        {
            auto const& entry = it->second;
            std::visit([&](auto* field) { action(*field, entry.name, entry.ownerKey); },
                       entry.get(model));
            return true;
        }
        return false;
    }

    std::unordered_map<std::string, Entry> fields;
};


/** \brief {name, computer} registry for MHD quantities that require an actual computation pass
 * ("/mhd/V", "/mhd/P") before they can be read via MHDAccessors - state.V/state.P are only kept
 * up to date at t=0 (see MHDState::initialize), so they need recomputing from the conservative
 * variables (rho, rhoV, B, Etot) on every dump. Mirrors FluidComputers; unlike it, the pressure
 * EOS needs gamma, which lives on the diagnostic (DiagnosticProperties::fileAttributes), not the
 * model, so compute() takes it as an explicit argument rather than pulling it off the state.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_mhd_model_v<Model>
struct ModelView<Hierarchy, Model>::MHDComputers
{
    using Computer = std::function<void(ModelView&, double, std::size_t, std::size_t, double)>;

    static MHDComputers& getOrCreateFor(auto& modelView)
    {
        if (!modelView.mhdComputers_)
            modelView.mhdComputers_ = std::make_unique<MHDComputers>(modelView);
        return *modelView.mhdComputers_;
    }

    MHDComputers(auto& modelView)
    {
        computers["/mhd/V"]
            = [](ModelView& mv, double /*gamma*/, std::size_t minLvl, std::size_t maxLvl,
                 double timestamp) { compute_mhd_velocity(mv, minLvl, maxLvl, timestamp); };
        computers["/mhd/P"]
            = [](ModelView& mv, double gamma, std::size_t minLvl, std::size_t maxLvl,
                 double timestamp) { compute_mhd_pressure(mv, gamma, minLvl, maxLvl, timestamp); };
    }

    // returns false if qty has no registered computer
    bool compute(std::string const& qty, ModelView& modelView, double gamma, std::size_t minLvl,
                 std::size_t maxLvl, double timestamp) const
    {
        auto it = computers.find(qty);
        if (it == computers.end())
            return false;
        it->second(modelView, gamma, minLvl, maxLvl, timestamp);
        return true;
    }

    std::unordered_map<std::string, Computer> computers;
};



} // namespace PHARE::diagnostic



#endif // DIAGNOSTIC_MODEL_VIEW_HPP
