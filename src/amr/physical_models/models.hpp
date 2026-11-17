#ifndef PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
#define PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP

#include "amr/resources_manager/amr_utils.hpp"

#include <tuple>
#include <string>
#include <cstddef>
#include <utility>
#include <type_traits>

namespace PHARE::amr
{

// patches of a level with the given resources set on each patch in turn
template<typename GridLayout, typename ResourcesManager, typename Level, typename... Args>
struct LevelAccessor
{
    using GridLayout_t = GridLayout;
    using level_t      = Level;

    LevelAccessor(ResourcesManager& rm, level_t& level, Args&... args_in)
        : rm_{rm}
        , level_{level}
        , args{args_in...}
    {
    }

    auto size() const { return static_cast<std::size_t>(level_.getLocalNumberOfPatches()); }

    // fn(view) per patch, the resources are set on the patch for the duration of the call
    void exec(auto&& fn)
    {
        std::apply(
            [&](auto&... resources) {
                for (auto const& patch : rm_.enumerate(level_, resources...))
                {
                    PatchView view{rm_, *patch, args};
                    fn(view);
                }
            },
            args);
    }

private:
    using patch_t = std::remove_cvref_t<decltype(*std::declval<level_t&>().getPatch(0))>;

    // only valid during exec, the patch is owned by the level
    struct PatchView
    {
        PatchView(ResourcesManager& rm, patch_t const& patch, std::tuple<Args&...> const& args_in)
            : id{to_string(patch.getGlobalId())}
            , layout{amr::layoutFromPatch<GridLayout_t>(patch)}
            , args{args_in}
            , rm_{rm}
            , patch_{patch}
        {
        }

        auto& patchID() const { return id; }

        void setTime(double const time, auto&... quantities) const
        {
            (rm_.setTime(quantities, patch_, time), ...);
        }

        // cells of this patch owned by another patch of the level, see makeNonOwnedBoxesFor
        auto nonOwnedBoxes(auto const& neighbours) const
        {
            return makeNonOwnedBoxesFor<GridLayout_t::dimension>(patch_, neighbours);
        }

        // this patch domain and its particle ghost box overlaps with the patches of the level, see
        // makeNonLevelGhostBoxFor
        auto nonLevelGhostBoxes(auto const& neighbours) const
        {
            return makeNonLevelGhostBoxFor<GridLayout_t>(patch_, neighbours);
        }

        // patch data not managed as a resource, e.g. the load balancer workload
        template<typename PatchData>
        auto& getPatchDataAs(int const id) const
        {
            return dynamic_cast<PatchData&>(*patch_.getPatchData(id));
        }

        std::string id;
        GridLayout_t layout;
        std::tuple<Args&...> args;

    private:
        ResourcesManager& rm_;
        patch_t const& patch_;
    };

public:
    using PatchView_t = PatchView;

private:
    ResourcesManager& rm_;
    level_t& level_;
    std::tuple<Args&...> args;
};

template<typename GridLayout, typename ResourcesManager, typename Level>
auto make_level_accessor(Level& level, ResourcesManager& rm, auto&... args)
{
    return LevelAccessor<GridLayout, ResourcesManager, Level,
                         std::remove_reference_t<decltype(args)>...>{rm, level, args...};
}

template<typename Model, typename Level>
auto make_model_level_accessor(Level& level, Model& model, auto&... args)
    requires std::is_same_v<std::remove_const_t<Level>, typename Model::amr_types::level_t>
{
    return make_level_accessor<typename Model::gridlayout_type>(level, *model.resourcesManager,
                                                                args...);
}

// fn(view) has run for every patch of the level on return
template<typename Accessor>
void onLevel(Accessor&& accessor, auto&& fn)
{
    accessor.exec(fn);
}

// sum of fn(view) over the patches of the level on this rank
template<typename Accessor>
auto sumOnLevel(Accessor&& accessor, auto&& fn)
{
    using PatchView_t = std::remove_cvref_t<Accessor>::PatchView_t;
    std::decay_t<decltype(fn(std::declval<PatchView_t&>()))> sum{};
    accessor.exec([&](auto& view) { sum += fn(view); });
    return sum;
}

} // namespace PHARE::amr

#endif // PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
