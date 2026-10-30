#ifndef DIAGNOSTIC_MODEL_VIEW_HPP
#define DIAGNOSTIC_MODEL_VIEW_HPP

#include "core/def.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/amr_constants.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/messengers/field_operate_transaction.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include "diagnostic/computers/diagnostic_computers.hpp"

#include "dict.hpp"

#include <SAMRAI/xfer/RefineAlgorithm.h>

#include <functional>
#include <memory>
#include <string>
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


template<typename Derived, typename Hierarchy, typename Model>
class BaseModelView : public IModelView
{
public:
    using GridLayout        = Model::gridlayout_type;
    using VecField          = Model::vecfield_type;
    using ResMan            = Model::resources_manager_type;
    using Field             = Model::field_type;
    using TensorFieldData_t = ResMan::template UserTensorField_t</*rank=*/2>::patch_data_type;
    static constexpr auto dimension = Model::dimension;

    using PatchProperties
        = cppdict::Dict<float, double, std::size_t, std::vector<int>, std::vector<std::uint32_t>,
                        std::vector<double>, std::vector<std::size_t>, std::string,
                        std::vector<std::string>>;

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
                                        std::forward<Action>(action), minLevel, maxLevel, *this,
                                        model_);
    }

    NO_DISCARD auto periodicities() const
    {
        auto const& periodicities = hierarchy_.periodicities();
        return std::vector<int>(periodicities.begin(), periodicities.end());
    }
    NO_DISCARD auto domainBox() const { return hierarchy_.domainBox(); }
    NO_DISCARD auto origin() const { return std::vector<double>(dimension, 0); }
    NO_DISCARD auto cellWidth() const { return hierarchy_.cellWidth(); }
    NO_DISCARD auto maxLevel() const { return hierarchy_.maxLevel(); }

    NO_DISCARD std::string getLayoutTypeString() const
    {
        return std::string{GridLayout::implT::type};
    }

    NO_DISCARD auto getPatchProperties(std::string patchID, GridLayout const& grid) const
    {
        PatchProperties dict;
        dict["origin"]   = grid.origin().toVector();
        dict["nbrCells"] = core::Point<std::uint32_t, Model::dimension>{grid.nbrCells()}.toVector();
        dict["lower"]    = grid.AMRBox().lower.toVector();
        dict["upper"]    = grid.AMRBox().upper.toVector();
        dict["mpi_rank"] = static_cast<std::size_t>(mpi::rank());
        return dict;
    }

    NO_DISCARD static auto getEmptyPatchProperties(PatchProperties dict = {})
    {
        dict["origin"]   = std::vector<double>{};
        dict["nbrCells"] = std::vector<std::uint32_t>{};
        dict["lower"]    = std::vector<int>{};
        dict["upper"]    = std::vector<int>{};
        dict["mpi_rank"] = std::size_t{0};
        return dict;
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


    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return derived().getCompileTimeResourcesViewList();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return derived().getCompileTimeResourcesViewList();
    }

protected:
    Model& model_;
    Hierarchy& hierarchy_;

private:
    Derived& derived() { return static_cast<Derived&>(*this); }
    Derived const& derived() const { return static_cast<Derived const&>(*this); }
};


template<typename Hierarchy, typename Model>
class ModelView;


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
class ModelView<Hierarchy, Model>
    : public BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>
{
    using Super        = BaseModelView<ModelView<Hierarchy, Model>, Hierarchy, Model>;
    using Field        = Model::field_type;
    using VecField     = Model::vecfield_type;
    using TensorFieldT = Model::ions_type::tensorfield_type;

public:
    using Model_t                = Model;
    using physical_quantity_type = Model::physical_quantity_type;

    struct FluidAccessors;
    struct FluidComputers;

    ModelView(Hierarchy& hierarchy, Model& model)
        : Super{hierarchy, model}
    {
        declareMomentumTensorAlgos();
    }

    NO_DISCARD VecField& getB() const { return this->model_.state.electromag.B; }

    NO_DISCARD VecField& getE() const { return this->model_.state.electromag.E; }

    NO_DISCARD auto& getIons() const { return this->model_.state.ions; }

    auto& tmpField() { return tmpField_; }

    auto& tmpVecField() { return tmpVec_; }

    template<std::size_t rank = 2>
    auto& tmpTensorField()
    {
        static_assert(rank > 0 and rank < 3);
        if constexpr (rank == 1)
            return tmpVec_;
        else
            return tmpTensor_;
    }

    void fillPopMomTensor(auto& lvl, auto const time, auto const popidx)
    {
        using value_type = TensorFieldT::value_type;
        auto constexpr N = core::detail::tensor_field_dim_from_rank<2>();

        auto& rm   = *this->model_.resourcesManager;
        auto& ions = this->model_.state.ions;

        for (auto patch : rm.enumerate(lvl, ions, tmpTensor_))
            for (std::uint8_t c = 0; c < N; ++c)
                std::memcpy(tmpTensor_[c].data(), ions[popidx].momentumTensor()[c].data(),
                            ions[popidx].momentumTensor()[c].size() * sizeof(value_type));

        MTAlgos[popidx].getOrCreateSchedule(this->hierarchy_, lvl.getLevelNumber()).fillData(time);

        for (auto patch : rm.enumerate(lvl, ions, tmpTensor_))
            for (std::uint8_t c = 0; c < N; ++c)
                std::memcpy(ions[popidx].momentumTensor()[c].data(), tmpTensor_[c].data(),
                            ions[popidx].momentumTensor()[c].size() * sizeof(value_type));
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(tmpField_, tmpVec_, tmpTensor_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(tmpField_, tmpVec_, tmpTensor_);
    }

protected:
    void declareMomentumTensorAlgos()
    {
        auto& rm = *this->model_.resourcesManager;

        auto const dst_name = tmpTensor_.name();

        for (auto& pop : this->model_.state.ions)
        {
            auto& MTAlgo        = MTAlgos.emplace_back();
            auto const src_name = pop.momentumTensor().name();

            auto&& [idDst, idSrc] = rm.getIDsList(dst_name, src_name);
            MTAlgo.MTalgo->registerRefine(
                idDst, idSrc, idDst, nullptr,
                std::make_shared<
                    amr::TensorFieldGhostInterpOverlapFillPattern<typename Super::GridLayout,
                                                                  /*rank_=*/2>>());
        }

        // can't create schedules here as the hierarchy has no levels yet
    }

    struct MTAlgo
    {
        auto& getOrCreateSchedule(auto& hierarchy, int const ilvl)
        {
            using PlusEqualsOp = core::PlusEquals<typename VecField::value_type>;

            auto const level   = hierarchy.getPatchLevel(ilvl);
            auto schedule_iter = MTschedules.find(ilvl);
            auto const create_schedule
                = schedule_iter == MTschedules.end() or schedule_iter->second.level.lock() != level;

            if (create_schedule)
                schedule_iter
                    = MTschedules
                          .insert_or_assign(
                              ilvl,
                              Entry{level,
                                    MTalgo->createSchedule(
                                        level, 0,
                                        std::make_shared<amr::FieldBorderOpTransactionFactory<
                                            typename Super::TensorFieldData_t, PlusEqualsOp>>())})
                          .first;

            return *schedule_iter->second.schedule;
        }

        struct Entry
        {
            std::weak_ptr<SAMRAI::hier::PatchLevel> level; // invalidadted if the Level is destroyed
            std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule;
        };

        std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> MTalgo
            = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();
        std::map<int, Entry> MTschedules;
    };

    std::vector<MTAlgo> MTAlgos;
    std::unique_ptr<FluidAccessors> fluidAccessors_;
    std::unique_ptr<FluidComputers> fluidComputers_;
    Field tmpField_{"PHARE_sumField", core::HybridQuantity::Scalar::rho};
    VecField tmpVec_{"PHARE_sumVec", core::HybridQuantity::Vector::V};
    TensorFieldT tmpTensor_{"PHARE_sumTensor", core::HybridQuantity::Tensor::M};
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

    // diag only
    NO_DISCARD VecField& getV() { return V_diag_; }

    NO_DISCARD const VecField& getV() const { return V_diag_; }

    NO_DISCARD Field& getP() { return P_diag_; }

    NO_DISCARD const Field& getP() const { return P_diag_; }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(V_diag_, P_diag_, tmpField_, tmpVec_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(V_diag_, P_diag_, tmpField_, tmpVec_);
    }

    auto& tmpField() { return tmpField_; }

    auto& tmpVecField() { return tmpVec_; }

    template<std::size_t rank = 2>
    auto& tmpTensorField()
    {
        static_assert(rank == 1);
        return tmpVec_;
    }

protected:
    // these quantities are not always up to date in the calculations but we can compute them from
    // the conservative variables when needed their registration and allocation are handled in the
    // model
    VecField V_diag_{"diagnostics_V_", core::MHDQuantity::Vector::V};
    Field P_diag_{"diagnostics_P_", core::MHDQuantity::Scalar::P};

    Field tmpField_{"PHARE_sumField_MHD", core::MHDQuantity::Scalar::ScalarAllPrimal};
    VecField tmpVec_{"PHARE_sumVec_MHD", core::MHDQuantity::Vector::VecAllPrimal};

    std::unique_ptr<MHDAccessors> mhdAccessors_;
    std::unique_ptr<MHDComputers> mhdComputers_;
};



/** \brief {name, accessor} registry for hybrid fluid quantities ("/ions/charge_density",
 * "/ions/pop/<name>/momentum_tensor", ...), keyed by the full diagnostic.quantity path. Built
 * once per ModelView from the real population names/indices, so callers do a single map lookup
 * instead of looping populations to find which one's tree prefix matches - shared across writer
 * types (h5, vtk) since it only depends on the ModelView, not on any writer.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
struct ModelView<Hierarchy, Model>::FluidAccessors
{
    using Field_t       = Model_t::field_type;
    using VecField_t    = Model_t::vecfield_type;
    using TensorField_t = Model_t::ions_type::tensorfield_type;
    using FieldVariant  = std::variant<Field_t*, VecField_t*, TensorField_t*>;
    using FieldFn       = std::function<FieldVariant(ModelView&)>;

    struct Entry
    {
        std::string ownerKey; // patchAttributes bucket: "ion" or "fluid_<pop name>"
        std::string name;     // bare quantity name, e.g. "density", "momentum_tensor"
        FieldFn get;
    };

    static FluidAccessors& getOrCreateFor(ModelView& modelView)
    {
        if (!modelView.fluidAccessors_)
            modelView.fluidAccessors_ = std::make_unique<FluidAccessors>(modelView);
        return *modelView.fluidAccessors_;
    }

    FluidAccessors(ModelView& modelView);

    // action(field, name, ownerKey) - false if qty has no registered accessor
    template<typename Action>
    bool dispatch(std::string const& qty, ModelView& modelView, Action&& action) const
    {
        if (auto it = fields.find(qty); it != fields.end())
        {
            auto const& entry = it->second;
            std::visit([&](auto* field) { action(*field, entry.name, entry.ownerKey); },
                       entry.get(modelView));
            return true;
        }
        return false;
    }

    // whether this model supports qty at all
    NO_DISCARD bool provides(std::string const& qty) const { return fields.count(qty) > 0; }

    std::unordered_map<std::string, Entry> fields;
};


template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
ModelView<Hierarchy, Model>::FluidAccessors::FluidAccessors(ModelView& modelView)
{
    auto& ions = modelView.getIons();

    fields["/ions/charge_density"]  = {"ion", "charge_density", [](ModelView& mv) -> FieldVariant {
                                          return &mv.getIons().chargeDensity();
                                      }};
    fields["/ions/mass_density"]    = {"ion", "mass_density", [](ModelView& mv) -> FieldVariant {
                                        return &mv.getIons().massDensity();
                                    }};
    fields["/ions/bulkVelocity"]    = {"ion", "bulkVelocity", [](ModelView& mv) -> FieldVariant {
                                        return &mv.getIons().velocity();
                                    }};
    fields["/ions/charge_flux"]     = {"ion", "charge_flux", [](ModelView& mv) -> FieldVariant {
                                       return &mv.getIons().chargeFlux();
                                   }};
    fields["/ions/momentum_tensor"] = {"ion", "momentum_tensor", [](ModelView& mv) -> FieldVariant {
                                           return &mv.getIons().momentumTensor();
                                       }};

    for (std::size_t idx = 0; idx < ions.size(); ++idx)
    {
        std::string const tree{"/ions/pop/" + ions[idx].name() + "/"};
        std::string const ownerKey{"fluid_" + ions[idx].name()};

        fields[tree + "density"] = {ownerKey, "density", [idx](ModelView& mv) -> FieldVariant {
                                        return &mv.getIons()[idx].particleDensity();
                                    }};
        fields[tree + "charge_density"]
            = {ownerKey, "charge_density",
               [idx](ModelView& mv) -> FieldVariant { return &mv.getIons()[idx].chargeDensity(); }};
        fields[tree + "flux"] = {ownerKey, "flux", [idx](ModelView& mv) -> FieldVariant {
                                     return &mv.getIons()[idx].flux();
                                 }};
        fields[tree + "momentum_tensor"]
            = {ownerKey, "momentum_tensor", [idx](ModelView& mv) -> FieldVariant {
                   return &mv.getIons()[idx].momentumTensor();
               }};
    }
}


/** \brief {name, computer} registry for hybrid fluid quantities that require an actual
 * computation pass ("/ions/momentum_tensor", "/ions/pop/<name>/momentum_tensor"), keyed by the
 * full diagnostic.quantity path. Writer-agnostic (takes minLvl/maxLvl/timestamp as plain
 * arguments rather than a concrete writer type) so it's shared across writer types the same way
 * FluidAccessors is.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_hybrid_model_v<Model>
struct ModelView<Hierarchy, Model>::FluidComputers
{
    using Computer = std::function<void(ModelView&, std::size_t, std::size_t, double)>;

    static FluidComputers& getOrCreateFor(ModelView& modelView)
    {
        if (!modelView.fluidComputers_)
            modelView.fluidComputers_ = std::make_unique<FluidComputers>(modelView);
        return *modelView.fluidComputers_;
    }

    FluidComputers(ModelView& modelView);

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
ModelView<Hierarchy, Model>::FluidComputers::FluidComputers(ModelView& modelView)
{
    computers["/ions/momentum_tensor"]
        = [](ModelView& mv, std::size_t minLvl, std::size_t maxLvl, double timestamp) {
              compute_momentum_tensor(mv, minLvl, maxLvl, timestamp);
          };

    auto& ions = modelView.getIons();
    for (std::size_t idx = 0; idx < ions.size(); ++idx)
        computers["/ions/pop/" + ions[idx].name() + "/momentum_tensor"]
            = [idx](ModelView& mv, std::size_t minLvl, std::size_t maxLvl, double timestamp) {
                  compute_pop_momentum_tensor(mv, idx, minLvl, maxLvl, timestamp);
              };
}


/** \brief {name, accessor} registry for MHD quantities ("/mhd/rho", "/mhd/V", ...), keyed by the
 * full diagnostic.quantity path - mirrors FluidAccessors, but the MHD state has no
 * per-population tree, so there is a single "mhd" bucket and one entry per field.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_mhd_model_v<Model>
struct ModelView<Hierarchy, Model>::MHDAccessors
{
    using Field_t      = Model_t::field_type;
    using VecField_t   = Model_t::vecfield_type;
    using FieldVariant = std::variant<Field_t*, VecField_t*>;
    using FieldFn      = std::function<FieldVariant(ModelView&)>;

    struct Entry
    {
        std::string ownerKey; // patchAttributes bucket - always "mhd"
        std::string name;     // bare quantity name, e.g. "rho", "V"
        FieldFn get;
    };

    static MHDAccessors& getOrCreateFor(ModelView& modelView)
    {
        if (!modelView.mhdAccessors_)
            modelView.mhdAccessors_ = std::make_unique<MHDAccessors>(modelView);
        return *modelView.mhdAccessors_;
    }

    // defined inline: clang's qualified lookup for an out-of-line nested ctor of this second
    // partial specialization fails to disambiguate an unqualified `ModelView` parameter type
    MHDAccessors(ModelView&)
    {
        std::string const tree{"/mhd/"};

        fields[tree + "rho"]
            = {"mhd", "rho", [](ModelView& mv) -> FieldVariant { return &mv.getRho(); }};
        fields[tree + "V"] = {"mhd", "V", [](ModelView& mv) -> FieldVariant { return &mv.getV(); }};
        fields[tree + "P"] = {"mhd", "P", [](ModelView& mv) -> FieldVariant { return &mv.getP(); }};
        fields[tree + "rhoV"]
            = {"mhd", "rhoV", [](ModelView& mv) -> FieldVariant { return &mv.getRhoV(); }};
        fields[tree + "Etot"]
            = {"mhd", "Etot", [](ModelView& mv) -> FieldVariant { return &mv.getEtot(); }};
    }

    // action(field, name, ownerKey) - false if qty has no registered accessor
    template<typename Action>
    bool dispatch(std::string const& qty, ModelView& modelView, Action&& action) const
    {
        if (auto it = fields.find(qty); it != fields.end())
        {
            auto const& entry = it->second;
            std::visit([&](auto* field) { action(*field, entry.name, entry.ownerKey); },
                       entry.get(modelView));
            return true;
        }
        return false;
    }

    // whether this model supports qty at all
    NO_DISCARD bool provides(std::string const& qty) const { return fields.count(qty) > 0; }

    std::unordered_map<std::string, Entry> fields;
};


/** \brief {name, computer} registry for MHD quantities that require an actual computation pass
 * ("/mhd/V", "/mhd/P") before they can be read via MHDAccessors - V/P are only kept up to date
 * at t=0, so they need recomputing from the conservative variables (rho, rhoV, B, Etot) on every
 * dump. Mirrors FluidComputers; unlike it, the pressure EOS needs gamma, which lives on the
 * diagnostic (DiagnosticProperties::fileAttributes), not the model, so compute() takes it as an
 * explicit argument.
 */
template<typename Hierarchy, typename Model>
    requires solver::is_mhd_model_v<Model>
struct ModelView<Hierarchy, Model>::MHDComputers
{
    using Computer = std::function<void(ModelView&, double, std::size_t, std::size_t, double)>;

    static MHDComputers& getOrCreateFor(ModelView& modelView)
    {
        if (!modelView.mhdComputers_)
            modelView.mhdComputers_ = std::make_unique<MHDComputers>(modelView);
        return *modelView.mhdComputers_;
    }

    MHDComputers(ModelView&)
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
