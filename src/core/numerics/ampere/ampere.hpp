#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include <core/utilities/types.hpp>

namespace PHARE::core
{

template<typename GridLayout>
class Ampere
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Ampere(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField& J)
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const shrink = ConstArray<std::uint32_t, dimension>(1);

        layout_.evalOnShrinkedGhostBox(
            Jx, shrink, [](auto&&... args) { JxEq_(args...); }, Jx, B, layout_);
        layout_.evalOnShrinkedGhostBox(
            Jy, shrink, [](auto&&... args) { JyEq_(args...); }, Jy, B, layout_);
        layout_.evalOnShrinkedGhostBox(
            Jz, shrink, [](auto&&... args) { JzEq_(args...); }, Jz, B, layout_);
    }


private:
    GridLayout layout_;

    static void JxEq_(auto const& ijk, auto& Jx, auto const& B, auto const& layout)
    {
        auto const& [_, By, Bz] = B();

        if constexpr (dimension == 1)
            Jx(ijk) = 0.0;

        if constexpr (dimension == 2)
            Jx(ijk) = layout.template deriv<Direction::Y>(Bz, ijk);

        if constexpr (dimension == 3)
            Jx(ijk) = layout.template deriv<Direction::Y>(Bz, ijk)
                      - layout.template deriv<Direction::Z>(By, ijk);
    }

    static void JyEq_(auto const& ijk, auto& Jy, auto const& B, auto const& layout)
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1 || dimension == 2)
            Jy(ijk) = -layout.template deriv<Direction::X>(Bz, ijk);

        if constexpr (dimension == 3)
            Jy(ijk) = layout.template deriv<Direction::Z>(Bx, ijk)
                      - layout.template deriv<Direction::X>(Bz, ijk);
    }

    static void JzEq_(auto const& ijk, auto& Jz, auto const& B, auto const& layout)
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1)
            Jz(ijk) = layout.template deriv<Direction::X>(By, ijk);

        else
            Jz(ijk) = layout.template deriv<Direction::X>(By, ijk)
                      - layout.template deriv<Direction::Y>(Bx, ijk);
    }
};

} // namespace PHARE::core
#endif
