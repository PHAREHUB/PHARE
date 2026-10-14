#ifndef PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
#define PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP

#include "core/utilities/thread_pool.hpp"
#include "core/data/field/field_tiles.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include <tuple>
#include <cstddef>
#include <utility>

namespace PHARE::amr
{

// Model is void when not made from a model (make_level_accessor)
template<typename GridLayout, typename Model, typename ResourcesManager, typename Level,
         typename... Args>
struct LevelAccessor
{
    using GridLayout_t = GridLayout;
    using Model_t      = Model;
    using level_t      = Level;

    LevelAccessor(ResourcesManager& rm, level_t& level, Args&... args_in)
        : rm_{rm}
        , level_{level}
        , n_patches{static_cast<std::size_t>(level_.getLocalNumberOfPatches())}
        , args{std::forward_as_tuple(args_in...)}
    {
    }

    auto size() const { return n_patches; }

    auto operator[](std::size_t const i)
    {
        return std::apply(
            [&](auto... resources) {
                auto looper = rm_.enumerate(level_, resources...);

                auto patch_view = looper[i];
                auto& patch     = *patch_view;
                return PatchView{core::to_string(patch.getGlobalId()),
                                 amr::layoutFromPatch<GridLayout_t>(patch),
                                 std::make_tuple(resources...)};
            },
            args);
    }

private:
    struct PatchView
    {
        auto& patchID() const { return id; }

        std::string id;
        GridLayout_t layout;
        std::tuple<Args...> args;
    };

    ResourcesManager& rm_;
    level_t& level_;
    std::size_t n_patches;
    std::tuple<Args&...> args;
};

template<typename GridLayout, typename Model = void, typename ResourcesManager, typename Level>
auto make_level_accessor(Level& level, ResourcesManager& rm, auto&... args)
{
    return LevelAccessor<GridLayout, Model, ResourcesManager, Level,
                         std::remove_reference_t<decltype(args)>...>{rm, level, args...};
}

template<typename Model>
auto make_model_level_accessor(typename Model::amr_types::level_t& level, Model& model,
                               auto&... args)
{
    return make_level_accessor<typename Model::gridlayout_type, Model>(
        level, *model.resourcesManager, args...);
}

// a pool per patch, a task per tile: fn(view, tile_idx) has run for every tile of every patch on
// return, tiled(view) gives the tile count
template<typename Accessor>
void tiled_level_exec(Accessor& accessor, auto const& tiled, auto const& fn)
{
    core::ThreadPool::INSTANCE().tile_exec(
        accessor.size(), [&](auto const i) { return accessor[i]; },
        [&](auto& view) { return core::tile_count(tiled(view)); }, fn);
}

// a pool per patch, one task each: fn(view) has run for every patch on return
template<typename Accessor>
void patch_level_exec(Accessor& accessor, auto const& fn)
{
    core::ThreadPool::INSTANCE().patch_exec(
        accessor.size(), [&](auto const i) { return accessor[i]; }, fn);
}

// fn(view, tile_idx) per tile if tiled(view) is tiled, else fn(view) per patch, so fn can use
// core::maybe_tile_at(x, tile_idx...) for either
template<typename Accessor>
void level_exec(Accessor& accessor, auto const& tiled, auto const& fn)
{
    using PatchView_t = decltype(accessor[0]);
    using Tiled_t     = std::remove_cvref_t<decltype(tiled(std::declval<PatchView_t&>()))>;

    if constexpr (core::is_field_tile_set_v<Tiled_t> or core::has_tiled_field_type_c<Tiled_t>)
        tiled_level_exec(accessor, tiled, fn);
    else
        patch_level_exec(accessor, fn);
}

} // namespace PHARE::amr

#endif // PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
