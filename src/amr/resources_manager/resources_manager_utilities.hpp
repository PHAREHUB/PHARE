#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_UTILITIES_HPP

#include "resources_meta.hpp"

#include <string>
#include <vector>
#include <type_traits>


namespace PHARE
{
namespace amr
{

    template<typename ResourceManager, typename ResourceView>
    class ResourceResolver
    {
        auto constexpr static resolve_t()
        {
            if constexpr (IsTensorFieldResource<ResourceView>)
                return typename ResourceManager::template UserTensorField_t<ResourceView::rank>{};
            else if constexpr (IsParticlesResource<ResourceView>)
                return typename ResourceManager::template UserParticle_t<ResourceView>{};
            else if constexpr (IsFieldResource<ResourceView>)
                return typename ResourceManager::UserField_t{};
            else
                throw std::runtime_error("bad condition");
        }

    public:
        using type = std::decay_t<decltype(resolve_t())>;

        auto static make_shared_variable(ResourceView const& view)
        {
            if constexpr (IsTensorFieldResource<ResourceView>)
                return std::make_shared<typename type::variable_type>(view.name(),
                                                                      view.physicalQuantity());
            else if constexpr (IsParticlesResource<ResourceView>)
                return std::make_shared<typename type::variable_type>(view.name());
            else if constexpr (IsFieldResource<ResourceView>)
                return std::make_shared<typename type::variable_type>(view.name(),
                                                                      view.physicalQuantity());
            else
                throw std::runtime_error("bad condition");
        }
    };



    /** extractNames of direct Field and Particle Resources of the given ResourceView
     * Is called by the other overload of extractNames()
     */
    template<typename ResourceView>
    void extractNames(ResourceView& view, std::vector<std::string>& names)
    {
        if constexpr (IsResource<ResourceView>)
            names.push_back(view.name());
    }


    /** @brief extractNames returns a vector of strings containing the names of all resources
     * associated with a ResourceView
     */
    template<typename ResourceView>
    std::vector<std::string> extractNames(ResourceView& view)
    {
        std::vector<std::string> names;

        if constexpr (HasCompileTimeSubResources<ResourceView>)
        {
            // unpack the tuple subResources and apply for each element registerResources()
            std::apply([&names](auto&... subResource) { (extractNames(subResource, names), ...); },
                       view.getCompileTimeResourceViewList());
        }

        if constexpr (HasRunTimeSubResources<ResourceView>)
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
