#ifndef PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
#define PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP


#include "amr/resources_manager/amr_utils.hpp"

#include <tuple>
#include <cstddef>


namespace PHARE::amr
{


template<typename Model, typename... Args>
struct ModelLevelAccessor
{
    using GridLayout_t = Model::gridlayout_type;
    using Model_t      = Model;
    using level_t      = Model::amr_types::level_t;

    ModelLevelAccessor(Model& model, level_t& level, Args&... args_in)
        : model_{model}
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
                auto looper = model_.resourcesManager->enumerate(level_, resources...);

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


    Model& model_;
    level_t& level_;
    std::size_t n_patches;
    std::tuple<Args...> args;
};



template<typename Model, typename... Args>
ModelLevelAccessor(Model&, typename Model::amr_types::level_t&, Args&...)
    -> ModelLevelAccessor<Model, Args...>;


template<typename Model>
auto make_model_level_accessor(typename Model::amr_types::level_t& level, Model& model,
                               auto&... args)
{
    return ModelLevelAccessor{model, level, args...};
}


} // namespace PHARE::amr


#endif // PHARE_AMR_PHYSICAL_MODEL_MODELS_HPP
