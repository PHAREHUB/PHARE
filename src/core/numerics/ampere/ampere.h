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
class Ampere : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField, typename... Boxes>
    void operator()(VecField const& B, VecField& J, Boxes&... boxes)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ampere - GridLayout not set, cannot proceed to calculate ampere()");

        auto const& [b0, b1, b2] = std::forward_as_tuple(boxes...);

        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        layout_->evalOnBox(b0, [&](auto&... args) mutable { JxEq_(Jx, B, args...); });
        layout_->evalOnBox(b1, [&](auto&... args) mutable { JyEq_(Jy, B, args...); });
        layout_->evalOnBox(b2, [&](auto&... args) mutable { JzEq_(Jz, B, args...); });
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField& J)
    {
        auto const& [Jx, Jy, Jz] = J();

        (*this)(B, J, Jx, Jy, Jz);
    }

    template<typename VecField, typename Box>
    static void op(GridLayout& layout, VecField const& B, VecField& J,
                   std::array<Box, 3> const& boxes)
    {
        Ampere self;
        self.setLayout(&layout);
        self(B, J, boxes[0], boxes[1], boxes[2]);
    }



private:
    template<typename VecField, typename Field, typename... Indexes>
    void JxEq_(Field& Jx, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [_, By, Bz] = B();

        if constexpr (dimension == 2)
            Jx(ijk...) = layout_->deriv(Bz, {ijk...}, DirectionTag<Direction::Y>{});

        if constexpr (dimension == 3)
            Jx(ijk...) = layout_->deriv(Bz, {ijk...}, DirectionTag<Direction::Y>{})
                         - layout_->deriv(By, {ijk...}, DirectionTag<Direction::Z>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void JyEq_(Field& Jy, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1)
            Jy(ijk...) = -layout_->deriv(Bz, {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 2)
            Jy(ijk...) = -layout_->deriv(Bz, {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 3)
            Jy(ijk...) = layout_->deriv(Bx, {ijk...}, DirectionTag<Direction::Z>{})
                         - layout_->deriv(Bz, {ijk...}, DirectionTag<Direction::X>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void JzEq_(Field& Jz, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [Bx, By, Bz] = B();

        if constexpr (dimension == 1)
            Jz(ijk...) = layout_->deriv(By, {ijk...}, DirectionTag<Direction::X>{});

        else
            Jz(ijk...) = layout_->deriv(By, {ijk...}, DirectionTag<Direction::X>{})
                         - layout_->deriv(Bx, {ijk...}, DirectionTag<Direction::Y>{});
    }
};

} // namespace PHARE::core
#endif
