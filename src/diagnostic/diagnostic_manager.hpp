#ifndef PHARE_DIAGNOSTIC_MANAGER_HPP_
#define PHARE_DIAGNOSTIC_MANAGER_HPP_

#include "amr/resources_manager/resources_manager_utilities.hpp"
#include "mpi/mpi_utils.hpp"

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"

#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/hybrid_model.hpp"

#include "initializer/data_provider.hpp"

#include "diagnostic_props.hpp"
#include "diagnostic_model_view.hpp"

#include <map>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <variant>
#include <utility>

namespace PHARE::diagnostic
{
enum class Mode { LIGHT, FULL };


class IDiagnosticsManager
{
public:
    virtual bool dump(double timeStamp)                          = 0;
    virtual bool dump(double timeStamp, double timeStep)         = 0;
    virtual void dump_level(std::size_t level, double timeStamp) = 0;
    virtual ~IDiagnosticsManager() {}
};


// finds whichever type among Models... is the hybrid (resp. mhd) model, regardless of
// position or how many models there are (one hybrid-only model, one mhd-only model, or
// both) - type is void if none match
template<typename... Models>
struct FindHybridModel
{
    using type = void;
};

template<typename Model, typename... Rest>
struct FindHybridModel<Model, Rest...>
{
    using type = std::conditional_t<solver::is_hybrid_model_v<Model>, Model,
                                    typename FindHybridModel<Rest...>::type>;
};

template<typename... Models>
struct FindMHDModel
{
    using type = void;
};

template<typename Model, typename... Rest>
struct FindMHDModel<Model, Rest...>
{
    using type = std::conditional_t<solver::is_mhd_model_v<Model>, Model,
                                    typename FindMHDModel<Rest...>::type>;
};

template<typename... Models>
struct HasModels
{
    using HybridModel_t = typename FindHybridModel<Models...>::type;
    using MHDModel_t    = typename FindMHDModel<Models...>::type;
};

template<typename Hierarchy, typename... Models>
struct DiagnosticsModelMapper
{
    static constexpr bool has_hybrid_model
        = !std::is_void_v<typename FindHybridModel<Models...>::type>;
    static constexpr bool has_mhd_model = !std::is_void_v<typename FindMHDModel<Models...>::type>;
    static constexpr auto dimension     = std::tuple_element_t<0, std::tuple<Models...>>::dimension;
    using This                          = DiagnosticsModelMapper<Hierarchy, Models...>;
    using HasModels_t                   = HasModels<Models...>;

    DiagnosticsModelMapper(auto& hier, auto&... models_)
        : hierarchy{hier}
        , models{std::forward_as_tuple(models_...)}
    {
        modelViews.reserve(sizeof...(models_)); // never reallocate - ModelView is never moved
        core::for_N<sizeof...(models_)>([&](auto i) {
            auto& model   = std::get<i>(this->models);
            using Model_t = std::decay_t<decltype(model)>;
            modelViews.emplace_back(std::in_place_type<ModelView<Hierarchy, Model_t>>, hier, model);
        });
    }

    NO_DISCARD auto cellWidth() const { return this->hierarchy.cellWidth(); }
    NO_DISCARD auto maxLevel() const { return this->hierarchy.maxLevel(); }

    // level enumeration is a hierarchy-level concern, not a per-model one - delegate straight
    // to the free function rather than resolving any particular model's ModelView
    template<typename... Args>
    void onLevels(Args&&... args)
    {
        amr::onLevels(this->hierarchy, std::forward<Args>(args)...);
    }


    using ModelViewVariant = std::variant<ModelView<Hierarchy, Models>...>;

    /** \brief visits [minLevel, maxLevel] across however many models are registered - hierarchy
     * knows which model owns each level (see Hierarchy::modelForLevel), so per level we pick the
     * matching ModelView and delegate to its own (single-model) visitHierarchy for that level.
     * ModelView itself stays single-model; only the writer needs to know about all of them.
     *
     * action is called as action(layout, patchID, levelNumber, modelViewVariant) - the variant
     * is the same one stored in modelViews (already holding the concrete, currently-patch-bound
     * ModelView for this callback), so callers that need model-specific accessors (getB(), etc.)
     * resolve it with std::visit/if constexpr rather than assuming a single fixed ModelView.
     */
    template<typename Action>
    void visitHierarchy(Action&& action, std::size_t minLevel, std::size_t maxLevel)
    {
        auto const hierLevels = static_cast<std::size_t>(this->hierarchy.getNumberOfLevels());

        for (std::size_t ilvl = minLevel; ilvl < hierLevels && ilvl <= maxLevel; ++ilvl)
        {
            auto const modelName = this->hierarchy.modelForLevel(ilvl);

            for (auto& modelViewVariant : modelViews)
                std::visit(
                    [&](auto& mv) {
                        using Model_t = typename std::decay_t<decltype(mv)>::Model_t;
                        if (Model_t::model_name == modelName)
                            mv.visitHierarchy(
                                [&](auto& layout, auto const& patchID, auto lvl) {
                                    action(layout, patchID, lvl, modelViewVariant);
                                },
                                ilvl, ilvl);
                    },
                    modelViewVariant);
        }
    }

    NO_DISCARD static constexpr bool is_hybrid_model(ModelViewVariant const& modelView)
        requires(has_hybrid_model)
    {
        return std::holds_alternative<ModelView<Hierarchy, typename HasModels_t::HybridModel_t>>(
            modelView);
    }

    auto& hyridModelView()
        requires(has_hybrid_model)
    {
        using HybridModelView_t = ModelView<Hierarchy, typename HasModels_t::HybridModel_t>;
        for (auto& mv : modelViews)
            if (std::holds_alternative<HybridModelView_t>(mv))
                return std::get<HybridModelView_t>(mv);
        throw std::runtime_error("hyridModelView: no hybrid model registered");
    }

    /** \brief resolves the ModelViewVariant that owns a given level, without visiting any
     * patches - for callers that need level (not patch) granularity, e.g. quantities valid
     * on every level regardless of model (electromag) where per-level init/write still needs
     * to know which concrete ModelView to dispatch to.
     */
    NO_DISCARD auto& levelModelView(std::size_t ilvl)
    {
        auto const modelName = this->hierarchy.modelForLevel(ilvl);
        for (auto& mv : modelViews)
        {
            bool const match = std::visit(
                [&](auto& concrete) {
                    using Model_t = typename std::decay_t<decltype(concrete)>::Model_t;
                    return Model_t::model_name == modelName;
                },
                mv);
            if (match)
                return mv;
        }
        throw std::runtime_error("levelModelView: no model registered for level "
                                 + std::to_string(ilvl));
    }


    NO_DISCARD static constexpr bool is_mhd_model(ModelViewVariant const& modelView)
        requires(has_mhd_model)
    {
        return std::holds_alternative<ModelView<Hierarchy, typename HasModels_t::MHDModel_t>>(
            modelView);
    }

    auto& mhdModelView()
        requires(has_mhd_model)
    {
        using MHDModelView_t = ModelView<Hierarchy, typename HasModels_t::MHDModel_t>;
        for (auto& mv : modelViews)
            if (std::holds_alternative<MHDModelView_t>(mv))
                return std::get<MHDModelView_t>(mv);
        throw std::runtime_error("mhdModelView: no mhd model registered");
    }

    template<typename Field_t>
        requires(amr::is_field_v<Field_t>)
    NO_DISCARD auto& tmpFor(Field_t const&)
    {
        using PQ = Field_t::physical_quantity_type;
        if constexpr (std::is_same_v<PQ, core::HybridQuantity::Scalar>)
            return hyridModelView().tmpField();
        else if constexpr (std::is_same_v<PQ, core::MHDQuantity::Scalar>)
            return mhdModelView().tmpField();
        else
            static_assert(core::dependent_false_v<Field_t>);
    }

    template<typename TensorField_t>
        requires(amr::is_tensor_field_v<TensorField_t>)
    NO_DISCARD auto& tmpFor(TensorField_t const&)
    {
        using PQ = TensorField_t::field_type::physical_quantity_type;
        if constexpr (std::is_same_v<PQ, core::HybridQuantity::Scalar>)
            return hyridModelView().template tmpTensorField<TensorField_t::rank>();
        else if constexpr (std::is_same_v<PQ, core::MHDQuantity::Scalar>)
            return mhdModelView().template tmpTensorField<TensorField_t::rank>();
        else
            static_assert(core::dependent_false_v<TensorField_t>);
    }

    Hierarchy& hierarchy;
    std::tuple<Models&...> models;
    std::vector<ModelViewVariant> modelViews;
};



template<typename Writer>
class DiagnosticsManager : public IDiagnosticsManager
{
public:
    // using Model_t = typename Writer::Model_t;

    bool dump(double timeStamp) override;
    bool dump(double timeStamp, double timeStep) override;


    void dump_level(std::size_t level, double timeStamp) override;


    // constructs writer_ in place from (hier, config, models...) - Writer is deliberately
    // non-movable (its H5TypeWriters hold a reference back to the owning Writer instance),
    // so it must never be built as a temporary and moved from
    DiagnosticsManager(auto& hier, typename Writer::Config_t const& config, auto&... models)
        : writer_{hier, config, models...}
    {
    }


    template<typename Hierarchy>
    NO_DISCARD static std::unique_ptr<DiagnosticsManager>
    make_unique(Hierarchy& hier, auto const& dict, auto&... models)
    {
        auto dMan
            = std::make_unique<DiagnosticsManager>(hier, Writer::Config_t::FROM(dict), models...);

        (dMan->registerDiagnostics(models, dict), ...);

        return dMan;
    }

    void registerDiagnostics(auto const& model, initializer::PHAREDict const& diagsParams);


    DiagnosticsManager& addDiagDict(initializer::PHAREDict const& dict);
    DiagnosticsManager& addDiagDict(initializer::PHAREDict&& dict) { return addDiagDict(dict); }


    NO_DISCARD auto& diagnostics() const { return diagnostics_; }


    NO_DISCARD auto& writer() { return writer_; }


    DiagnosticsManager(DiagnosticsManager const&)            = delete;
    DiagnosticsManager(DiagnosticsManager&&)                 = delete;
    DiagnosticsManager& operator=(DiagnosticsManager const&) = delete;
    DiagnosticsManager& operator=(DiagnosticsManager&&)      = delete;

private:
    NO_DISCARD bool needsAction_(double nextTime, double timeStamp, double timeStep)
    {
        // casting to float to truncate double to avoid trailing imprecision
        return static_cast<float>(std::abs(nextTime - timeStamp)) < static_cast<float>(timeStep);
    }

    bool needsElapsedAction_(double const nextTime) const
    {
        return mpi::unix_timestamp_now() > nextTime;
    }



    NO_DISCARD bool needsWrite_(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto const& diag_key     = diag.type + diag.quantity;
        auto& nextWriteTimestamp = nextWrite_[diag_key];
        auto& nextWriteElapsed   = nextWriteElapsed_[diag_key];

        auto const writeTimestampNow
            = nextWriteTimestamp < diag.writeTimestamps.size()
              and needsAction_(diag.writeTimestamps[nextWriteTimestamp], timeStamp, timeStep);

        auto const writeElapsedNow
            = nextWriteElapsed < diag.elapsedTimestamps.size()
              and needsElapsedAction_(diag.elapsedTimestamps[nextWriteElapsed]);

        if (writeTimestampNow)
            ++nextWriteTimestamp;
        if (writeElapsedNow)
            ++nextWriteElapsed;

        return writeTimestampNow || writeElapsedNow;
    }


    NO_DISCARD bool needsCompute_(DiagnosticProperties& diag, double timeStamp, double timeStep)
    {
        auto nextCompute = nextCompute_[diag.type + diag.quantity];
        return nextCompute < diag.computeTimestamps.size()
               and needsAction_(diag.computeTimestamps[nextCompute], timeStamp, timeStep);
    }


    Writer writer_;
    std::map<std::string, std::size_t> nextCompute_;
    std::map<std::string, std::size_t> nextWrite_;
    std::map<std::string, std::size_t> nextWriteElapsed_;
    std::vector<DiagnosticProperties> diagnostics_;

    std::time_t const start_time_{mpi::unix_timestamp_now()};
};

template<typename Model_t>
auto constexpr diagnostic_types_per_model(Model_t const& model)
{
    if constexpr (solver::is_hybrid_model_v<Model_t>)
        return std::vector<std::string>{"fluid", "electromag", "particle", "meta", "info"};

    else if constexpr (solver::is_mhd_model_v<Model_t>)
        return std::vector<std::string>{"mhd", "meta", "electromag"};

    else
        static_assert(core::dependent_false_v<Model_t>, "Unsupported model type");

    throw std::runtime_error("registerDiagnostics: Shouldn't happen!");
}


template<typename ModelMapper>
void DiagnosticsManager<ModelMapper>::registerDiagnostics(auto const& model,
                                                          initializer::PHAREDict const& diagsParams)
{
    for (auto& diagType : diagnostic_types_per_model(model))
    {
        // several diags of the same type can be registeroed
        // corresponding to several blocks in user input
        // fluid0, fluid1, electromag0, electromag1, etc.
        std::size_t diagBlockID = 0;
        while (diagsParams.contains(diagType)
               && diagsParams[diagType].contains(diagType + std::to_string(diagBlockID)))
        {
            std::string const diagName = diagType + std::to_string(diagBlockID);
            addDiagDict(diagsParams[diagType][diagName]);
            ++diagBlockID;
        }
    }
}


template<typename ModelMapper>
DiagnosticsManager<ModelMapper>&
DiagnosticsManager<ModelMapper>::addDiagDict(initializer::PHAREDict const& diagParams)
{
    auto& diagProps           = diagnostics_.emplace_back(DiagnosticProperties{});
    diagProps.type            = diagParams["type"].template to<std::string>();
    diagProps.quantity        = diagParams["quantity"].template to<std::string>();
    diagProps.writeTimestamps = diagParams["write_timestamps"].template to<std::vector<double>>();

    if (diagParams.contains("elapsed_timestamps"))
    {
        diagProps.elapsedTimestamps
            = diagParams["elapsed_timestamps"].template to<std::vector<double>>();
        for (auto& time : diagProps.elapsedTimestamps)
            time += start_time_; // expected for comparison later
    }

    diagProps["flush_every"] = diagParams["flush_every"].template to<std::size_t>();

    diagProps.computeTimestamps
        = diagParams["compute_timestamps"].template to<std::vector<double>>();

    diagProps.nAttributes = diagParams["n_attributes"].template to<std::size_t>();
    for (std::size_t i = 0; i < diagProps.nAttributes; ++i)
    {
        std::string const idx = std::to_string(i);
        std::string const key = diagParams["attribute_" + idx + "_key"];
        diagProps.forward_file_attribute(key, diagParams["attribute_" + idx + "_value"]);
    }

    return *this;
}


template<typename ModelMapper>
void DiagnosticsManager<ModelMapper>::dump_level(std::size_t level, double timeStamp)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;

    for (auto& diag : diagnostics_)
        activeDiagnostics.emplace_back(&diag);

    // compute() now happens inside writer_.dump()/dump_level() itself, immediately before
    // each diagnostic's own write pass - see H5Writer::writeDatasets_
    writer_.dump_level(level, activeDiagnostics, timeStamp);
}


template<typename ModelMapper>
bool DiagnosticsManager<ModelMapper>::dump(double timeStamp, double timeStep)
{
    std::vector<DiagnosticProperties*> activeDiagnostics;
    for (auto& diag : diagnostics_)
    {
        auto diagID = diag.type + diag.quantity;

        if (needsWrite_(diag, timeStamp, timeStep))
            activeDiagnostics.emplace_back(&diag);
    }

    if (activeDiagnostics.size() > 0)
    {
        PHARE_LOG_SCOPE(1, "DiagnosticsManager::dump(double, double)");
        writer_.dump(activeDiagnostics, timeStamp);
    }

    return activeDiagnostics.size() > 0;
}


template<typename ModelMapper>
bool DiagnosticsManager<ModelMapper>::dump(double const timeStamp)
{
    std::vector<DiagnosticProperties*> diagnostics;
    for (auto& diag : diagnostics_)
        diagnostics.emplace_back(&diag);

    if (diagnostics.size() > 0)
    {
        PHARE_LOG_SCOPE(1, "DiagnosticsManager::dump()");
        writer_.dump(diagnostics, timeStamp);
    }

    return diagnostics.size() > 0;
}

} // namespace PHARE::diagnostic

#endif /* PHARE_DIAGNOSTIC_MANAGER_HPP_ */
