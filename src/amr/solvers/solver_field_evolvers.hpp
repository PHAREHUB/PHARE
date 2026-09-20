#ifndef PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP

#include "core/utilities/thread_pool.hpp"
#include "core/data/field/field_tiles.hpp"
#include "core/numerics/ampere/ampere.hpp"
#include "core/numerics/faraday/faraday.hpp"

#include "amr/physical_models/models.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include <memory>
#include <vector>

namespace PHARE::solver
{

// a pool per patch, a task per tile: fn(view, tile_idx) has run for every tile of every patch on
// return, tiled(view) gives the tile count
template<typename PatchView_t>
void tiled_level_exec(std::vector<std::shared_ptr<PatchView_t>> const& views, auto const& tiled,
                      auto const& fn)
{
    auto& tp = core::ThreadPool::INSTANCE();
    for (auto const& view : views)
    {
        auto& pool         = tp.get_pool(tp.first_ready_idx());
        auto const n_tiles = core::tile_count(tiled(*view));
        for (std::size_t t = 0; t < n_tiles; ++t)
            pool.detach_task([=] { fn(*view, t); });
    }
    tp.sync();
}

template<typename Accessor>
auto make_patch_views(Accessor& accessor)
{
    using PatchView_t = decltype(accessor[0]);
    std::vector<std::shared_ptr<PatchView_t>> views;
    views.reserve(accessor.size());
    for (std::size_t i = 0; i < accessor.size(); ++i)
        views.emplace_back(std::make_shared<PatchView_t>(accessor[i]));
    return views;
}

template<typename Accessor>
void tiled_level_exec(Accessor& accessor, auto const& tiled, auto const& fn)
{
    tiled_level_exec(make_patch_views(accessor), tiled, fn);
}

// on_tile(view, tile_idx) runs for every tile of every patch before the ghosts of output(view)
// are synced, also per tile, as syncing a tile reads its neighbours' domains
template<typename Accessor>
void tiled_level_transform(Accessor& accessor, auto const& output, auto const& on_tile)
{
    auto const views = make_patch_views(accessor);
    tiled_level_exec(views, output, on_tile);
    tiled_level_exec(views, output, [=](auto& view, auto const tile_idx) {
        core::sync_inner_ghosts(output(view), tile_idx);
    });
}

class FaradaySingleTransformer
{
    template<typename GridLayout>
    void operate(GridLayout const& layout, auto&&... args)
    {
        core::Faraday<GridLayout>{layout}(args...);
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E, VecField& Bnew,
                    double dt)
        requires(not core::has_tiled_field_type_c<VecField>)
    {
        operate(layout, B, E, Bnew, dt);
    }

    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& /*layout*/, VecField const& B, VecField const& E,
                    VecField& Bnew, double dt)
        requires(core::has_tiled_field_type_c<VecField>)
    {
        for (std::size_t i = 0; i < core::tile_count(Bnew); ++i)
            on_tile(i, B, E, Bnew, dt);

        core::sync_inner_ghosts(Bnew);
    }

    // one tile's worth of work, Bnew ghosts are not synced
    void on_tile(std::size_t const tile_idx, auto const& B, auto const& E, auto& Bnew,
                 double const dt)
    {
        core::tile_exec_with_layout_at(
            tile_idx, [&](auto& layout, auto&&... args) { operate(layout, args...); }, B, E, Bnew,
            dt);
    }
};

template<typename Model>
class FaradayLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

public:
    explicit FaradayLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField const& E, VecField& Bnew,
                    double dt)
    {
        FaradaySingleTransformer{}(layout, B, E, Bnew, dt);
    }

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
                FaradaySingleTransformer{}.on_tile(tile_idx, B_v, E_v, Bnew_v, dt);
            });
    }

    level_t& level_;
    Model& model_;
};

template<typename Model>
FaradayLevelTransformer(typename Model::amr_types::level_t&, Model&)
    -> FaradayLevelTransformer<Model>;


class AmpereSingleTransformer
{
    template<typename GridLayout>
    void operate(GridLayout const& layout, auto&&... args)
    {
        core::Ampere<GridLayout>{layout}(args...);
    }

public:
    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J)
        requires(not core::has_tiled_field_type_c<VecField>)
    {
        operate(layout, B, J);
    }

    template<typename GridLayout, typename VecField>
    void operator()(GridLayout const& /*layout*/, VecField const& B, VecField& J)
        requires(core::has_tiled_field_type_c<VecField>)
    {
        for (std::size_t i = 0; i < core::tile_count(J); ++i)
            on_tile(i, B, J);

        core::sync_inner_ghosts(J);
    }

    // one tile's worth of work, J ghosts are not synced
    void on_tile(std::size_t const tile_idx, auto const& B, auto& J)
    {
        core::tile_exec_with_layout_at(
            tile_idx, [&](auto& layout, auto&&... args) { operate(layout, args...); }, B, J);
    }
};

template<typename Model>
class AmpereLevelTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

public:
    explicit AmpereLevelTransformer(level_t& level, auto& model)
        : level_{level}
        , model_{model}
    {
    }

    template<typename VecField>
    void operator()(GridLayout const& layout, VecField const& B, VecField& J)
    {
        AmpereSingleTransformer{}(layout, B, J);
    }

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
                AmpereSingleTransformer{}.on_tile(tile_idx, B_v, J_v);
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
    using Faraday_t = FaradayLevelTransformer<Model>;
    using Ampere_t  = AmpereLevelTransformer<Model>;
};

} // namespace PHARE::solver

#endif /* PHARE_AMR_SOLVERS_SOLVER_FIELD_EVOLVERS_HPP */
