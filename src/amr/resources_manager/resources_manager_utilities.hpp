#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP

#include "core/utilities/meta/meta_utilities.hpp"

#include <string>
#include <type_traits>
#include <vector>

namespace PHARE
{
namespace amr
{


    template<typename ResourcesUser, typename Attempt = void>
    struct has_particles : std::false_type
    {
    };

    template<typename ResourcesUser, typename Attempt = void>
    struct has_sub_resources : std::false_type
    {
    };


    template<typename ResourcesUser, typename Attempt = void>
    struct has_runtime_subresourceuser_list : std::false_type
    {
    };


    template<typename ResourcesUser, typename Attempt = void>
    struct has_compiletime_subresourcesuser_list : std::false_type
    {
    };



    /** \brief is_field is a traits that permit to check if a ResourcesUser
     * is a field
     */
    template<typename ResourcesUser, typename Attempt = void>
    struct is_field : std::false_type
    {
    };

    template<typename ResourcesUser>
    struct is_field<ResourcesUser, core::tryToInstanciate<
                                       decltype(std::declval<ResourcesUser>().physicalQuantity())>>
        : std::true_type
    {
    };



    /** \brief has_particles is a traits that permit to check if a ResourcesUser
     * has particles
     */
    template<typename ResourcesUser>
    struct has_particles<
        ResourcesUser,
        core::tryToInstanciate<decltype(std::declval<ResourcesUser>().getParticleArrayNames())>>
        : std::true_type
    {
    };



    /** @brief has_runtime_subresourceuser_list is a compile-time function that returns true if the
     * given ResourcesUser has a runtime list of ResourcesUsers, like a vector of ResourcesUsers.
     */
    template<typename ResourcesUser>
    struct has_runtime_subresourceuser_list<
        ResourcesUser, core::tryToInstanciate<
                           decltype(std::declval<ResourcesUser>().getRunTimeResourcesUserList())>>
        : std::true_type
    {
    };


    /** @brief has_compiletime_subresourcesuser_list is a compile-time function that returns true if
     * the given ResourcesUser has one or several ResourcesUsers that can be put in a compile-time
     * list.
     */
    template<typename ResourcesUser>
    struct has_compiletime_subresourcesuser_list<
        ResourcesUser, core::tryToInstanciate<decltype(std::declval<ResourcesUser>()
                                                           .getCompileTimeResourcesUserList())>>
        : std::true_type
    {
    };




    /** UseResourcePtr is used to select the resources patch data */
    struct UseResourcePtr
    {
    };


    /** UseNullPtr is used to select a nullptr with the correct type */
    struct UseNullPtr
    {
    };



    /** extractNames of direct Field and Particle Resources of the given ResourcesUser
     * Is called by the other overload of extractNames()
     */
    template<typename ResourcesUser>
    void extractNames(ResourcesUser& user, std::vector<std::string>& names)
    {
        if constexpr (is_field<ResourcesUser>::value)
        {
            names.push_back(user.getFieldNameAndQuantity().name);
        }

        if constexpr (has_particles<ResourcesUser>::value)
        {
            auto pnames = user.getParticleArrayNames();
            for (auto const& p : pnames)
            {
                names.push_back(p.name);
            }
        }
    }


    /** @brief extractNames returns a vector of strings containing the names of all resources
     * associated with a ResourcesUser
     */
    template<typename ResourcesUser>
    std::vector<std::string> extractNames(ResourcesUser& user)
    {
        std::vector<std::string> names;

        if constexpr (has_compiletime_subresourcesuser_list<ResourcesUser>::value)
        {
            // get a tuple here
            auto&& subResources = user.getCompileTimeResourcesUserList();

            // unpack the tuple subResources and apply for each element registerResources()
            std::apply([&names](auto&... subResource) { (extractNames(subResource, names), ...); },
                       subResources);
        }

        if constexpr (has_runtime_subresourceuser_list<ResourcesUser>::value)
        {
            auto&& resourcesUsers = user.getRunTimeResourcesUserList();
            for (auto& resourcesUser : resourcesUsers)
            {
                extractNames(resourcesUser, names);
            }
        }

        extractNames(user, names);

        return names;
    }


} // namespace amr

} // namespace PHARE

#endif // PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
