#ifndef PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_HPP
#define PHARE_CORE_DATA_VECFIELD_VECFIELD_COMPONENT_HPP

#include <unordered_map>
#include <string>
#include <stdexcept>
#include <cstdint>

namespace PHARE
{
namespace core
{
    enum class Component : std::uint16_t { X = 0, Y, Z, XX, XY, XZ, YY, YZ, ZZ };

    struct Components
    {
        template<std::size_t rank = 1>
        static std::unordered_map<std::string, Component> const& componentMap()
        {
            if constexpr (rank == 1)
                return vecComponentMap;
            else if constexpr (rank == 2)
                return tensComponentMap;
            else
                throw std::runtime_error("rank must be 1 or 2");
        }
        template<typename T>
        static Component at(T t)
        {
            return componentMap().at(std::forward<T>(t));
        }

        static std::unordered_map<std::string, Component> const vecComponentMap;
        static std::unordered_map<std::string, Component> const tensComponentMap;

        template<auto Tag>
        constexpr static bool check()
        {
            return Tag == Component::X || Tag == Component::Y || Tag == Component::Z
                   || Tag == Component::XX || Tag == Component::XY || Tag == Component::XZ
                   || Tag == Component::YY || Tag == Component::YZ || Tag == Component::ZZ;
        }
    };
    inline std::unordered_map<std::string, Component> const Components::vecComponentMap{
        {"x", Component::X}, {"y", Component::Y}, {"z", Component::Z}};
    inline std::unordered_map<std::string, Component> const Components::tensComponentMap{
        {"xx", Component::XX}, {"xy", Component::XY}, {"xz", Component::XZ},
        {"yy", Component::YY}, {"yz", Component::YZ}, {"zz", Component::ZZ}};

    struct VectorComponents
    {
        auto static map() { return Components::componentMap(); }
    };

    struct TensorComponents
    {
        auto static map() { return Components::componentMap<2>(); }
    };

} // namespace core
} // namespace PHARE

#endif
