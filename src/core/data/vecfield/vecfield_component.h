
#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_H
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_H

#include <unordered_map>
#include <string>

namespace PHARE
{
namespace core
{
    enum class Component { X, Y, Z };

    template<Component value>
    struct ComponentTag
    {
    };

    template<>
    struct ComponentTag<Component::X>
    {
        static const auto component = Component::X;
    };

    template<>
    struct ComponentTag<Component::Y>
    {
        static const auto component = Component::Y;
    };

    template<>
    struct ComponentTag<Component::Z>
    {
        static const auto component = Component::Z;
    };

    struct Components
    {
        template<typename T>
        static Component at(T t)
        {
            return componentMap.at(std::forward<T>(t));
        }

        static std::unordered_map<std::string, Component> const componentMap;
    };
    inline std::unordered_map<std::string, Component> const Components::componentMap{
        {"x", Component::X}, {"y", Component::Y}, {"z", Component::Z}};

} // namespace core
} // namespace PHARE

#endif
