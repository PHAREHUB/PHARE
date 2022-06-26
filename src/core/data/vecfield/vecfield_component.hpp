
#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_HPP
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_HPP

#include <unordered_map>
#include <string>

namespace PHARE
{
namespace core
{
    enum class Component { X, Y, Z };


    struct Components
    {
        template<typename T>
        static Component at(T t)
        {
            return componentMap.at(std::forward<T>(t));
        }

        static std::unordered_map<std::string, Component> const componentMap;

        template<auto Tag>
        constexpr static bool check()
        {
            return Tag == Component::X || Tag == Component::Y || Tag == Component::Z;
        }
    };
    inline std::unordered_map<std::string, Component> const Components::componentMap{
        {"x", Component::X}, {"y", Component::Y}, {"z", Component::Z}};

} // namespace core
} // namespace PHARE

#endif
