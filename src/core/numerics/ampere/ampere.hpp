#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include <core/utilities/types.hpp>

#include <cstdint>

namespace PHARE::core
{


enum class AmpereMode { ShrinkedGhost, GrownPhysical };

struct AmpereBox // structural type — usable as NTTP
{
    AmpereMode mode      = AmpereMode::ShrinkedGhost;
    std::uint32_t amount = 1;
};


template<typename GridLayout, AmpereBox box = AmpereBox{}>
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

        auto const amount = ConstArray<std::size_t, dimension>(box.amount);

        auto const eval = [&](auto& Jc, auto&& fn) {
            if constexpr (box.mode == AmpereMode::ShrinkedGhost)
                layout_.evalOnShrinkedGhostBox(Jc, amount, fn, Jc, B, layout_);
            else
                layout_.evalOnBiggerBox(Jc, amount, fn, Jc, B, layout_);
        };

        eval(Jx, [](auto&&... args) { JxEq_(args...); });
        eval(Jy, [](auto&&... args) { JyEq_(args...); });
        eval(Jz, [](auto&&... args) { JzEq_(args...); });
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
