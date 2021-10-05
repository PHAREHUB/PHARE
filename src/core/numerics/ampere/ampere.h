#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_H
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_H

#include <cstddef>
#include <iostream>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/utilities/index/index.h"


namespace PHARE::core
{
template<typename GridLayout>
struct StandardAmpereComputer
{
    constexpr static auto dimension = GridLayout::dimension;

    template<typename VecField, typename Field, typename... Indexes>
    void Jx(Field& Jx, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [_, By, Bz] = B();

        if constexpr (dimension == 2)
            Jx(ijk...) = layout.deriv(Bz, {ijk...}, DirectionTag<Direction::Y>{});

        if constexpr (dimension == 3)
            Jx(ijk...) = layout.deriv(Bz, {ijk...}, DirectionTag<Direction::Y>{})
                         - layout.deriv(By, {ijk...}, DirectionTag<Direction::Z>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void Jy(Field& Jy, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1)
            Jy(ijk...) = -layout.deriv(Bz, {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 2)
            Jy(ijk...) = -layout.deriv(Bz, {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 3)
            Jy(ijk...) = layout.deriv(Bx, {ijk...}, DirectionTag<Direction::Z>{})
                         - layout.deriv(Bz, {ijk...}, DirectionTag<Direction::X>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void Jz(Field& Jz, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1)
            Jz(ijk...) = layout.deriv(By, {ijk...}, DirectionTag<Direction::X>{});

        else
            Jz(ijk...) = layout.deriv(By, {ijk...}, DirectionTag<Direction::X>{})
                         - layout.deriv(Bx, {ijk...}, DirectionTag<Direction::Y>{});
    }

    GridLayout& layout;
};

template<typename GridLayout, typename Computer = StandardAmpereComputer<GridLayout>>
class Ampere : public LayoutHolder<GridLayout>
{
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField>
    void operator()(VecField const& B_, VecField& J_)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ampere - GridLayout not set, cannot proceed to calculate ampere()");

        auto B = B_.as_view();
        auto J = J_.as_view();

        Computer op{*this->layout_};

        layout_->scan(J_(Component::X),
                      [=](auto const&... args) mutable { op.Jx(J(Component::X), B, args...); });
        layout_->scan(J_(Component::Y),
                      [=](auto const&... args) mutable { op.Jy(J(Component::Y), B, args...); });
        layout_->scan(J_(Component::Z),
                      [=](auto const&... args) mutable { op.Jz(J(Component::Z), B, args...); });
    }
};

} // namespace PHARE::core
#endif
