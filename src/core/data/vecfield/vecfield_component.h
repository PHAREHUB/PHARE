
#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_H
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_H

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

} // namespace core
} // namespace PHARE

#endif
