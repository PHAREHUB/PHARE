#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include "core/data/ions/ion_population/particle_pack.hpp"

#include "field_resource.hpp"
#include "particle_resource.hpp"


#include <string>
#include <vector>
#include <type_traits>


namespace PHARE
{
namespace amr
{

    /** \brief is_field is a traits that permit to check if a ResourceView
     * is a field
     */
    template<typename ResourceView, typename Attempt = void>
    struct is_field : std::false_type
    {
    };

    template<typename ResourcesUser>
    struct is_field<ResourcesUser, core::tryToInstanciate<
                                       decltype(std::declval<ResourcesUser>().physicalQuantity())>>
        : std::true_type
    {
    };
    template<typename ResourceView>
    bool constexpr static is_field_v = is_field<ResourceView>::value;


    /** \brief is_tensor_field is a trait to check if a ResourceView is a tensor field
     */
    template<typename ResourceView, typename Attempt = void>
    struct is_tensor_field : std::false_type
    {
    };

    template<typename ResourcesUser>
    struct is_tensor_field<
        ResourcesUser, core::tryToInstanciate<decltype(std::declval<ResourcesUser>().components())>>
        : std::true_type
    {
    };
    template<typename ResourceView>
    bool constexpr static is_tensor_field_v = is_tensor_field<ResourceView>::value;


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

    template<typename ResourceManager, typename ResourceView>
    class ResourceResolver
    {
        auto constexpr static resolve_t()
        {
            if constexpr (is_tensor_field_v<ResourceView>)
                return typename ResourceManager::template UserTensorField_t<ResourceView::rank>{};
            else if constexpr (is_particles_v<ResourceView>)
                return typename ResourceManager::template UserParticle_t<ResourceView>{};
            else if constexpr (is_field_v<ResourceView>)
                return typename ResourceManager::UserField_t{};
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
