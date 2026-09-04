#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP

#include "core/core_meta.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include "core/data/ions/ion_population/particle_pack.hpp"

#include "field_resource.hpp"
#include "particle_resource.hpp"
#include "tensor_field_resource.hpp"


#include <string>
#include <tuple>
#include <vector>
#include <type_traits>


namespace PHARE
{
namespace amr
{

    using core::is_field_v;
    using core::is_tensor_field_v;


    /** \brief is_particles is a traits that permit to check if a ResourceView
     * has particles
     */
    template<typename ResourceView, typename Attempt = void>
    struct is_particles : std::false_type
    {
    };

    template<typename ResourceView>
    struct is_particles<
        ResourceView,
        core::tryToInstanciate<decltype(std::declval<ResourceView>().setBuffer(
            static_cast<core::ParticlesPack<typename ResourceView::particle_array_type>*>(
                nullptr)))>> : std::true_type
    {
    };
    template<typename ResourceView>
    bool constexpr static is_particles_v = is_particles<ResourceView>::value;



    template<typename ResourceView>
    struct is_resource
    {
        bool constexpr static value
            = core::any(is_field_v<ResourceView>, is_tensor_field_v<ResourceView>,
                        is_particles_v<ResourceView>);
    };
    template<typename ResourceView>
    bool constexpr static is_resource_v = is_resource<ResourceView>::value;



    /** \brief ResourcesManagerTypes bundles, for a single model's core_types (e.g.
     * core::PHARE_Types<opts>::Hybrid / ::MHD), the SAMRAI variable/patchdata types used to
     * register that model's resources.
     */
    template<typename core_types>
    struct ResourcesManagerTypes
    {
        using Grid_t       = core_types::Grid_t;
        using GridLayoutT  = core_types::GridLayout_t;
        using QuantityType = decltype(GridLayoutT::options.field_options)::Quantity;

        using UserField_t = UserFieldType<Grid_t, GridLayoutT>;

        template<std::size_t rank>
        using UserTensorField_t = UserTensorFieldType<rank, Grid_t, GridLayoutT, QuantityType>;

        template<typename ResourcesView>
        using UserParticle_t = UserParticleType<ResourcesView, GridLayoutT>;
    };


    /** \brief tells whether core_types is the entry (out of a ResourcesManager's TypeTuple) that
     * belongs to the given Model (HybridModel/MHDModel, ...) - matched on GridLayout_t/Grid_t,
     * which are unique per model.
     */
    template<typename Model, typename core_types>
    constexpr bool is_model_core_types()
    {
        return std::is_same_v<typename Model::gridlayout_type, typename core_types::GridLayout_t>
               and std::is_same_v<typename Model::grid_type, typename core_types::Grid_t>;
    }

    template<typename Model, typename... core_types_per_model>
    struct find_model_core_types
    {
        using type = void;
    };

    template<typename Model, typename core_types0, typename... core_types_rest>
    struct find_model_core_types<Model, core_types0, core_types_rest...>
    {
        using type = std::conditional_t<
            is_model_core_types<Model, core_types0>(), core_types0,
            typename find_model_core_types<Model, core_types_rest...>::type>;
    };

    template<typename Model, typename TypeTuple>
    struct model_core_types_of;

    template<typename Model, typename... core_types_per_model>
    struct model_core_types_of<Model, std::tuple<core_types_per_model...>>
    {
        using type = typename find_model_core_types<Model, core_types_per_model...>::type;
    };

    /** \brief resolves the ResourcesManagerTypes (UserField_t/UserTensorField_t/UserParticle_t)
     * that belong to Model, by looking Model up in TypeTuple (a ResourcesManager::TypeTuple) -
     * i.e. resolving *through* the (possibly multi-model) ResourcesManager rather than
     * independently re-deriving Model's shape.
     */
    template<typename Model, typename TypeTuple>
    struct ResourceUserTypesFor_
    {
        using core_types = typename model_core_types_of<Model, TypeTuple>::type;
        static_assert(!std::is_void_v<core_types>,
                      "Model is not one of the models registered on this ResourcesManager");
        using type = ResourcesManagerTypes<core_types>;
    };


    /** \brief has_ions is a trait that tells whether a model's core_types exposes particles
     * (i.e. an Ions_t) - only some models (e.g. Hybrid) do.
     */
    template<typename core_types>
    concept HasIons = requires { typename core_types::Ions_t; };

    template<typename core_types>
    bool constexpr static has_ions_v = HasIons<core_types>;


    /** \brief owns_resource tells whether the given model's core_types is the one that owns the
     * given (already concrete, model-specific) ResourceView.
     *
     * Fields and TensorFields carry their own physical quantity value (e.g. HybridQuantity::Scalar
     * or MHDQuantity::Vector), which is an exact, unambiguous match against exactly one model's
     * QuantityType::Scalar/Vector/Tensor. Particles carry no such tag, so models without an
     * Ions_t are simply excluded, and the remaining match is done on particle_array_type.
     */
    template<typename core_types, typename ResourceView>
    constexpr bool owns_resource()
    {
        if constexpr (is_field_v<ResourceView> or is_tensor_field_v<ResourceView>)
        {
            using Quantity = typename ResourcesManagerTypes<core_types>::QuantityType;
            using Qty_t    = std::decay_t<decltype(std::declval<ResourceView>().physicalQuantity())>;

            return std::is_same_v<Qty_t, typename Quantity::Scalar>
                   or std::is_same_v<Qty_t, typename Quantity::Vector>
                   or std::is_same_v<Qty_t, typename Quantity::Tensor>;
        }
        else if constexpr (is_particles_v<ResourceView>)
        {
            if constexpr (has_ions_v<core_types>)
                return std::is_same_v<typename ResourceView::particle_array_type,
                                      typename core_types::Ions_t::particle_array_type>;
            else
                return false;
        }
        else
            return false;
    }


    template<typename ResourceView, typename... core_types_per_model>
    struct find_owner_model
    {
        using type = void;
    };

    template<typename ResourceView, typename core_types0, typename... core_types_rest>
    struct find_owner_model<ResourceView, core_types0, core_types_rest...>
    {
        using type = std::conditional_t<
            owns_resource<core_types0, ResourceView>(), core_types0,
            typename find_owner_model<ResourceView, core_types_rest...>::type>;
    };

    template<typename ResourceView, typename TypeTuple>
    struct owner_model_of;

    template<typename ResourceView, typename... core_types_per_model>
    struct owner_model_of<ResourceView, std::tuple<core_types_per_model...>>
    {
        using type = typename find_owner_model<ResourceView, core_types_per_model...>::type;
    };


    template<typename ResourceManager, typename ResourceView>
    class ResourceResolver
    {
        using OwnerModel_t
            = typename owner_model_of<ResourceView, typename ResourceManager::TypeTuple>::type;

        static_assert(!std::is_void_v<OwnerModel_t>,
                      "ResourceView does not belong to any model registered on this "
                      "ResourcesManager");

        using RMT_t = ResourcesManagerTypes<OwnerModel_t>;

        auto constexpr static resolve_t()
        {
            if constexpr (is_tensor_field_v<ResourceView>)
                return typename RMT_t::template UserTensorField_t<ResourceView::rank>{};
            else if constexpr (is_particles_v<ResourceView>)
                return typename RMT_t::template UserParticle_t<ResourceView>{};
            else if constexpr (is_field_v<ResourceView>)
                return typename RMT_t::UserField_t{};
            else
                throw std::runtime_error("bad condition");
        }

    public:
        using type = std::decay_t<decltype(resolve_t())>;

        auto static make_shared_variable(ResourceView const& view)
        {
            if constexpr (is_tensor_field_v<ResourceView>)
                return std::make_shared<typename type::variable_type>(view.name(),
                                                                      view.physicalQuantity());
            else if constexpr (is_particles_v<ResourceView>)
                return std::make_shared<typename type::variable_type>(view.name());
            else if constexpr (is_field_v<ResourceView>)
                return std::make_shared<typename type::variable_type>(view.name(),
                                                                      view.physicalQuantity());
            else
                throw std::runtime_error("bad condition");
        }
    };


    /** @brief has_runtime_subresourceview_list is a compile-time function that returns true if the
     * given ResourceView has a runtime list of ResourceViews, like a vector of ResourceViews.
     */
    template<typename ResourceView, typename Attempt = void>
    struct has_runtime_subresourceview_list : std::false_type
    {
    };

    template<typename ResourceView>
    struct has_runtime_subresourceview_list<
        ResourceView, core::tryToInstanciate<
                          decltype(std::declval<ResourceView>().getRunTimeResourcesViewList())>>
        : std::true_type
    {
    };
    template<typename ResourceView>
    bool constexpr static has_runtime_subresourceview_list_v
        = has_runtime_subresourceview_list<ResourceView>::value;


    /** @brief has_compiletime_subresourcesview_list is a compile-time function that returns true if
     * the given ResourceView has one or several ResourceViews that can be put in a compile-time
     * list.
     */
    template<typename ResourceView, typename Attempt = void>
    struct has_compiletime_subresourcesview_list : std::false_type
    {
    };
    template<typename ResourceView>
    struct has_compiletime_subresourcesview_list<
        ResourceView, core::tryToInstanciate<
                          decltype(std::declval<ResourceView>().getCompileTimeResourcesViewList())>>
        : std::true_type
    {
    };
    template<typename ResourceView>
    bool constexpr static has_compiletime_subresourcesview_list_v
        = has_compiletime_subresourcesview_list<ResourceView>::value;


    template<typename RV>
    struct has_sub_resources
    {
        bool constexpr static value
            = has_compiletime_subresourcesview_list_v<RV> or has_runtime_subresourceview_list_v<RV>;
    };
    template<typename ResourceView>
    bool constexpr static has_sub_resources_v = has_sub_resources<ResourceView>::value;



    /** extractNames of direct Field and Particle Resources of the given ResourceView
     * Is called by the other overload of extractNames()
     */
    template<typename ResourceView>
    void extractNames(ResourceView& view, std::vector<std::string>& names)
    {
        if constexpr (is_resource_v<ResourceView>)
            names.push_back(view.name());
    }


    /** @brief extractNames returns a vector of strings containing the names of all resources
     * associated with a ResourceView
     */
    template<typename ResourceView>
    std::vector<std::string> extractNames(ResourceView& view)
    {
        std::vector<std::string> names;

        if constexpr (has_compiletime_subresourcesview_list<ResourceView>::value)
        {
            // unpack the tuple subResources and apply for each element registerResources()
            std::apply([&names](auto&... subResource) { (extractNames(subResource, names), ...); },
                       view.getCompileTimeResourceViewList());
        }

        if constexpr (has_runtime_subresourceview_list<ResourceView>::value)
        {
            for (auto& resourcesUser : view.getRunTimeResourceViewList())
                extractNames(resourcesUser, names);
        }

        extractNames(view, names);

        return names;
    }


} // namespace amr

} // namespace PHARE

#endif // PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
