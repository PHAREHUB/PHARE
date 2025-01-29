#ifndef PHARE_AMR_TOOLS_RESOURCES_GUARDS_HPP
#define PHARE_AMR_TOOLS_RESOURCES_GUARDS_HPP

#include "core/def/phare_mpi.hpp"

#include "resources_manager_utilities.hpp"

#include <memory>
#include <tuple>
#include <utility>

#include <SAMRAI/hier/Patch.h>

namespace PHARE
{
namespace amr
{
    /** \brief ViewGuard ... TODO fix doc
     */
    template<typename ResourcesManager, typename... Views>
    class ViewGuard
    {
    public:
        ViewGuard(SAMRAI::hier::Patch const& patch, ResourcesManager const& resourcesManager,
                  Views&... views)
            : views_{views...}
            , patch_{patch}
            , resourcesManager_{resourcesManager}
        {
            std::apply(
                [this](auto&... user) {
                    ((resourcesManager_.setResources_(user, UseResourcePtr{}, patch_)), ...);
                },
                views_);
        }




        ~ViewGuard()
        {
            std::apply(
                [this](auto&... view) {
                    ((resourcesManager_.setResources_(view, UseNullPtr{}, patch_)), ...);
                },
                views_);
        }




        ViewGuard(ViewGuard&&)                        = default;
        ViewGuard()                                   = delete;
        ViewGuard(ViewGuard const&)                   = delete;
        ViewGuard& operator=(ViewGuard const& source) = delete;
        ViewGuard& operator=(ViewGuard&&)             = delete;



    private:
        std::tuple<Views&...> views_;
        SAMRAI::hier::Patch const& patch_;
        ResourcesManager const& resourcesManager_;
    };
} // namespace amr
} // namespace PHARE
#endif
