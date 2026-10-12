#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

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

        Point<std::uint32_t, dimension> amount;

        for (size_t i = 0; i < dimension; ++i)
        {
            amount[i] = box.amount;
        }

        auto eval = [&](auto& Jc, auto&& fn) {
            if constexpr (box.mode == AmpereMode::ShrinkedGhost)
                layout_.evalOnShrinkedGhostBox(Jc, amount, fn);
            else
                layout_.evalOnBiggerBox(Jc, amount, fn);
        };

        eval(Jx, [&](auto&... args) mutable { JxEq_(Jx, B, args...); });
        eval(Jy, [&](auto&... args) mutable { JyEq_(Jy, B, args...); });
        eval(Jz, [&](auto&... args) mutable { JzEq_(Jz, B, args...); });
    }


private:
    GridLayout layout_;


    template<typename VecField, typename Field, typename... Indexes>
    void JxEq_(Field& Jx, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [_, By, Bz] = B();

        if constexpr (dimension == 1)
            Jx(ijk...) = 0.0;

        if constexpr (dimension == 2)
            Jx(ijk...) = layout_.template deriv<Direction::Y>(Bz, {ijk...});

        if constexpr (dimension == 3)
            Jx(ijk...) = layout_.template deriv<Direction::Y>(Bz, {ijk...})
                         - layout_.template deriv<Direction::Z>(By, {ijk...});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void JyEq_(Field& Jy, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1 || dimension == 2)
            Jy(ijk...) = -layout_.template deriv<Direction::X>(Bz, {ijk...});

        if constexpr (dimension == 3)
            Jy(ijk...) = layout_.template deriv<Direction::Z>(Bx, {ijk...})
                         - layout_.template deriv<Direction::X>(Bz, {ijk...});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void JzEq_(Field& Jz, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1)
            Jz(ijk...) = layout_.template deriv<Direction::X>(By, {ijk...});

        else
            Jz(ijk...) = layout_.template deriv<Direction::X>(By, {ijk...})
                         - layout_.template deriv<Direction::Y>(Bx, {ijk...});
    }
};

} // namespace PHARE::core
#endif
