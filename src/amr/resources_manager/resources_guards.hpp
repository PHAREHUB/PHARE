#ifndef PHARE_AMR_TOOLS_RESOURCES_GUARDS_HPP
#define PHARE_AMR_TOOLS_RESOURCES_GUARDS_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "resources_manager_utilities.hpp"

#include <memory>
#include <tuple>
#include <utility>

#include <SAMRAI/hier/Patch.h>

namespace PHARE
{
namespace amr
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
            std::apply(
                [this](auto&... user) { ((resourcesManager_.setResources_(user, patch_)), ...); },
                resourcesUsers_);
        }




        ~ResourcesGuard()
        {
            // set nullptr to all users in resourcesUsers_
            std::apply([this](auto&... user) { ((resourcesManager_.unsetResources_(user)), ...); },
                       resourcesUsers_);
        }




        // We just need the move constructor for using ResourceManager::createResourcesGuards
        ResourcesGuard(ResourcesGuard&&)                        = default;
        ResourcesGuard()                                        = delete;
        ResourcesGuard(ResourcesGuard const&)                   = delete;
        ResourcesGuard& operator=(ResourcesGuard const& source) = delete;
        ResourcesGuard& operator=(ResourcesGuard&&)             = delete;



    private:
        std::tuple<ResourcesUsers&...> resourcesUsers_;
        SAMRAI::hier::Patch const& patch_;
        ResourcesManager const& resourcesManager_;
    };
} // namespace amr
} // namespace PHARE
#endif
