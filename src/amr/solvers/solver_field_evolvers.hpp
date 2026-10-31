#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/data/field/field_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/physical_models/models.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/physical_models/hybrid_model.hpp"


namespace PHARE::solver
{

using amr::tiled_level_exec;

// on_tile(view, tile_idx) runs for every tile of every patch before the ghosts of output(view)
// are synced, also per tile, as syncing a tile reads its neighbours' domains
template<typename Accessor>
void tiled_level_transform(Accessor& accessor, auto const& output, auto const& on_tile)
{
    tiled_level_exec(accessor, output, on_tile);
    tiled_level_exec(accessor, output, [&](auto& view, auto const tile_idx) {
        core::sync_inner_ghosts(output(view), tile_idx);
    });
}


template<typename Model>
class FaradayLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Faraday<GridLayout>;

public:
    explicit FaradayLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    template<typename VecField>
    void operator()(VecField& B, VecField& E, VecField& Bnew, double const dt)
        requires(not core::has_tiled_field_type_c<VecField>)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, B, E, Bnew))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, E, Bnew, dt);
        }
    }

    template<typename VecField>
    void operator()(VecField& B, VecField& E, VecField& Bnew, double const dt)
        requires(core::has_tiled_field_type_c<VecField>)
    {
        auto accessor = amr::make_model_level_accessor(level_, model_, B, E, Bnew);
        tiled_level_transform(
            accessor, [](auto& view) -> auto& { return std::get<2>(view.args); },
            [dt](auto& view, auto const tile_idx) {
                auto& [B_v, E_v, Bnew_v] = view.args;
                core::tile_exec_with_layout_at(
                    tile_idx,
                    [](auto& layout, auto&&... args) {
                        using Layout_t = std::decay_t<decltype(layout)>;
                        core::Faraday<Layout_t>{layout}(args...);
                    },
                    B_v, E_v, Bnew_v, dt);
            });
    }

    level_t& level_;
    Model& model_;
};
template<typename Model>
FaradayLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> FaradayLevelTransformer<Model>;




template<typename Model, core::AmpereBox box = core::AmpereBox{}>
class AmpereLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::Ampere<GridLayout, box>;

public:
    explicit AmpereLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    void operator()(GridLayout& layout, auto&&... args) { core_type{layout}(args...); }

    template<typename VecField>
    void operator()(VecField& B, VecField& J)
        requires(not core::has_tiled_field_type_c<VecField>)
    {
        auto& rm = *model_.resourcesManager;
        for (auto& patch : rm.enumerate(level_, B, J))
        {
            auto layout = amr::layoutFromPatch<GridLayout>(*patch);
            (*this)(layout, B, J);
        }
    }

    template<typename VecField>
    void operator()(VecField& B, VecField& J)
        requires(core::has_tiled_field_type_c<VecField>)
    {
        auto accessor = amr::make_model_level_accessor(level_, model_, B, J);
        tiled_level_transform(
            accessor, [](auto& view) -> auto& { return std::get<1>(view.args); },
            [](auto& view, auto const tile_idx) {
                auto& [B_v, J_v] = view.args;
                core::tile_exec_with_layout_at(
                    tile_idx,
                    [](auto& layout, auto&&... args) {
                        using Layout_t = std::decay_t<decltype(layout)>;
                        core::Ampere<Layout_t, box>{layout}(args...);
                    },
                    B_v, J_v);
            });
    }

    level_t& level_;
    Model& model_;
};


template<typename Model>
AmpereLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> AmpereLevelTransformer<Model>;




template<typename level_t, typename Model>
struct TimeSetter
{
    void operator()(auto&... quantities)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, quantities...))
            (model.resourcesManager->setTime(quantities, *patch, newTime), ...);
    }

    level_t& level;
    Model& model;
    double newTime;
};

template<typename level_t, typename Model>
TimeSetter(level_t&, Model&, double) -> TimeSetter<level_t, Model>;



template<typename Model>
struct FieldEvolverDispatchers
{
    // hybrid: J is never communicated — Ampere computes it on physical+1,
    // exactly the layer Ohm reads. MHD keeps ghostbox-1 (its working state).
    static constexpr core::AmpereBox ampere_box
        = is_hybrid_model_v<Model> ? core::AmpereBox{core::AmpereMode::GrownPhysical, 1}
                                   : core::AmpereBox{};

    using Faraday_t = FaradayLevelTransformer<Model>;
    using Ampere_t  = AmpereLevelTransformer<Model, ampere_box>;
};


} // namespace PHARE::solver



#endif /* PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP */
