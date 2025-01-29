#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP

#include "core/utilities/meta/meta_utilities.hpp"

#include "field_resource.hpp"
#include "particle_resource.hpp"
#include "core/data/ions/ion_population/particle_pack.hpp"

#include <string>
#include <type_traits>
#include <vector>


namespace PHARE
{
namespace amr
{

    /** \brief is_field is a traits that permit to check if a View
     * is a field
     */
    template<typename View, typename Attempt = void>
    struct is_field : std::false_type
    {
    };

    template<typename View>
    struct is_field<View, core::tryToInstanciate<decltype(std::declval<View>().physicalQuantity())>>
        : std::true_type
    {
    };
    template<typename View>
    bool constexpr static is_field_v = is_field<View>::value;


    /** \brief is_particles is a traits that permit to check if a View
     * has particles
     */
    template<typename View, typename Attempt = void>
    struct is_particles : std::false_type
    {
    };

    template<typename View>
    struct is_particles<
        View, core::tryToInstanciate<decltype(std::declval<View>().setBuffer(
                  static_cast<core::ParticlesPack<typename View::particle_array_type>*>(nullptr)))>>
        : std::true_type
    {
    };
    template<typename View>
    bool constexpr static is_particles_v = is_particles<View>::value;



    template<typename View>
    struct is_final_view
    {
        bool constexpr static value = is_field_v<View> or is_particles_v<View>;
    };
    template<typename View>
    bool constexpr static is_final_view_v = is_final_view<View>::value;


    template<typename ResourceManager, typename View>
    class ViewInfoResolver
    {
        auto constexpr static resolve_t()
        {
            if constexpr (is_field_v<View>)
                return typename ResourceManager::Super::FieldViewInfo_t{};

            else if constexpr (is_particles_v<View>)
                return typename ResourceManager::Super::ParticleViewInfo_t{};

            else
                throw std::runtime_error("bad condition");
        }

    public:
        using view_info_type = std::decay_t<decltype(resolve_t())>;

        auto static make_shared_variable(View const& view)
        {
            if constexpr (is_field_v<View>)
                return std::make_shared<typename view_info_type::variable_type>(
                    view.name(), view.physicalQuantity());
            else
                return std::make_shared<typename view_info_type::variable_type>(view.name());
        }
    };


    /** @brief has_runtime_subresourceview_list is a compile-time function that returns true if the
     * given View has a runtime list of Views, like a vector of Views.
     */
    template<typename View, typename Attempt = void>
    struct has_runtime_views : std::false_type
    {
    };

    template<typename View>
    struct has_runtime_views<
        View, core::tryToInstanciate<decltype(std::declval<View>().getRunTimeResourcesViewList())>>
        : std::true_type
    {
    };
    template<typename View>
    bool constexpr static has_runtime_views_v = has_runtime_views<View>::value;


    /** @brief has_compiletime_subresourcesview_list is a compile-time function that returns true if
     * the given View has one or several Views that can be put in a compile-time
     * list.
     */
    template<typename View, typename Attempt = void>
    struct has_compiletime_views : std::false_type
    {
    };
    template<typename View>
    struct has_compiletime_views<
        View,
        core::tryToInstanciate<decltype(std::declval<View>().getCompileTimeResourcesViewList())>>
        : std::true_type
    {
    };
    template<typename View>
    bool constexpr static has_compiletime_views_v = has_compiletime_views<View>::value;


    template<typename View>
    struct has_subview
    {
        bool constexpr static value = has_compiletime_views_v<View> or has_runtime_views_v<View>;
    };
    template<typename View>
    bool constexpr static has_subview_v = has_subview<View>::value;

    /** UseResourcePtr is used to select the resources patch data */
    struct UseResourcePtr
    {
    };


    /** UseNullPtr is used to select a nullptr with the correct type */
    struct UseNullPtr
    {
    };



    /** extractNames of direct Field and Particle Resources of the given View
     * Is called by the other overload of extractNames()
     */
    template<typename View>
    void extractNames(View& view, std::vector<std::string>& names)
    {
        if constexpr (is_final_view_v<View>)
            names.push_back(view.name());
    }


    /** @brief extractNames returns a vector of strings containing the names of all resources
     * associated with a View
     */
    template<typename View>
    std::vector<std::string> extractNames(View& view)
    {
        std::vector<std::string> names;

        if constexpr (has_compiletime_views_v<View>)
        {
            // unpack the tuple subResources and apply for each element registerResources()
            std::apply([&names](auto&... views) { (extractNames(views, names), ...); },
                       view.getCompileTimeViewList());
        }

        if constexpr (has_runtime_views_v<View>)
        {
            for (auto& rt_view : view.getRunTimeViewList())
                extractNames(rt_view, names);
        }

        extractNames(view, names);

        return names;
    }


} // namespace amr

} // namespace PHARE

#endif // PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
