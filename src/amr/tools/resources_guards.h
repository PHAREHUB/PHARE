#ifndef PHARE_AMR_TOOLS_RESOURCES_GUARDS_H
#define PHARE_AMR_TOOLS_RESOURCES_GUARDS_H

#include "resources_manager_utilities.h"

#include <memory>
#include <tuple>
#include <utility>

#include <SAMRAI/hier/Patch.h>

namespace PHARE
{
/** \brief ResourcesGuards maintain a link with several resources in order to give them
 * an access to the data.
 *
 *  At construction it will take the patch, and all resourcesManagerUser object that need to
 *  be set via the ResourcesManager. Upon destruction, it will put all the previous
 *  object in an inactive state (for now it just put nullptr on them)

 * TODO: add the active mechanism
 */
template<typename ResourcesManager, typename... ResourcesUsers>
class ResourcesGuard
{
private:
    template<typename RequestedPtr, typename ResourcesUser>
    void setResources_(RequestedPtr requestedPtr, ResourcesUser& resourcesUser) const
    {
        resourcesManager_.setResources_(resourcesUser, requestedPtr, patch_);
    }
    template<typename RequestedPtr, typename ResourcesUser, typename... ResourcesUserTail>
    void setResources_(RequestedPtr requestedPtr, ResourcesUser& resourcesUser,
                       ResourcesUserTail&... resourcesUsers) const
    {
        setResources_(requestedPtr, resourcesUser);
        setResources_(requestedPtr, resourcesUsers...);
    }


    template<std::size_t... I>
    void resetResources_(std::index_sequence<I...>) const
    {
        setResources_(UseNullPtr{}, std::get<I>(resourcesUsers_)...);
    }


public:
    /**
     *  \brief At construction, each resourcesUser will be valid. At destruction time, each
     * resourcesUser will be invalid
     *
     *  \param[in] patch
     *  \param[in] resourcesManager
     *  \param[in,out] resourcesUsers
     */
    ResourcesGuard(SAMRAI::hier::Patch const& patch, ResourcesManager const& resourcesManager,
                   ResourcesUsers&... resourcesUsers)
        : resourcesUsers_{resourcesUsers...}
        , patch_{patch}
        , resourcesManager_{resourcesManager}
    {
        setResources_(UseResourcePtr{}, resourcesUsers...);
    }
    ~ResourcesGuard()
    {
        // Here we want to do a loop on the tuple : resourcesUsers_
        // we will do it in a recursive way, first with indexes we get a vector like that
        // 0 , 1 , 2 , ...
        // then we call setResources with the variadic arguments obtained by unpacking the tuple
        // with std::get
        auto indexes
            = std::make_index_sequence<std::tuple_size<decltype(resourcesUsers_)>::value>{};
        resetResources_(indexes);
    }

    // We just need the move constructor for using ResourceManager::createResourcesGuards
    ResourcesGuard(ResourcesGuard&&) = default;

    ResourcesGuard()                      = delete;
    ResourcesGuard(ResourcesGuard const&) = delete;
    ResourcesGuard& operator=(ResourcesGuard const& source) = delete;
    ResourcesGuard& operator=(ResourcesGuard&&) = delete;

private:
    std::tuple<ResourcesUsers&...> resourcesUsers_;
    SAMRAI::hier::Patch const& patch_;
    ResourcesManager const& resourcesManager_;
};
} // namespace PHARE
#endif
