#ifndef PHARE_FARADAY_H
#define PHARE_FARADAY_H

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/grid/gridlayout_utils.h"
#include "core/data/vecfield/vecfield_component.h"


namespace PHARE::core
{
template<typename GridLayout>
class Faraday : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField, typename... Boxes>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt,
                    Boxes&... boxes)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Faraday - GridLayout not set, cannot proceed to calculate faraday()");

        if (!(B.isUsable() && E.isUsable() && Bnew.isUsable()))
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        this->dt_ = dt;

        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        auto const& [b0, b1, b2] = std::forward_as_tuple(boxes...);

        layout_->evalOnBox(b0, [&](auto&... args) mutable { BxEq_(Bx, E, Bxnew, args...); });
        layout_->evalOnBox(b1, [&](auto&... args) mutable { ByEq_(By, E, Bynew, args...); });
        layout_->evalOnBox(b2, [&](auto&... args) mutable { BzEq_(Bz, E, Bznew, args...); });
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
    {
        if (!Bnew.isUsable())
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        (*this)(B, E, Bnew, dt, Bxnew, Bynew, Bznew);
    }

    template<typename VecField, typename Box>
    static void op(GridLayout& layout, VecField const& B, VecField const& E, VecField& Bnew,
                   double dt, std::array<Box, 3> const& boxes)
    {
        Faraday self;
        self.setLayout(&layout);
        self(B, E, Bnew, dt, boxes[0], boxes[1], boxes[2]);
    }


private:
    double dt_;


    template<typename VecField, typename Field, typename... Indexes>
    void BxEq_(Field const& Bx, VecField const& E, Field& Bxnew, Indexes const&... ijk) const
    {
        auto const& [_, Ey, Ez] = E();

        if constexpr (dimension == 1)
            Bxnew(ijk...) = Bx(ijk...);

        if constexpr (dimension == 2)
            Bxnew(ijk...)
                = Bx(ijk...)
                  - dt_ * layout_->deriv(E(Component::Z), {ijk...}, DirectionTag<Direction::Y>{});

        if constexpr (dimension == 3)
            Bxnew(ijk...) = Bx(ijk...)
                            - dt_ * layout_->deriv(Ez, {ijk...}, DirectionTag<Direction::Y>{})
                            + dt_ * layout_->deriv(Ey, {ijk...}, DirectionTag<Direction::Z>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void ByEq_(Field const& By, VecField const& E, Field& Bynew, Indexes const&... ijk) const
    {
        auto const& [Ex, _, Ez] = E();

        if constexpr (dimension == 1)
            Bynew(ijk...)
                = By(ijk...)
                  + dt_ * layout_->deriv(E(Component::Z), {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 2)
            Bynew(ijk...)
                = By(ijk...)
                  + dt_ * layout_->deriv(E(Component::Z), {ijk...}, DirectionTag<Direction::X>{});

        if constexpr (dimension == 3)
            Bynew(ijk...) = By(ijk...)
                            - dt_ * layout_->deriv(Ex, {ijk...}, DirectionTag<Direction::Z>{})
                            + dt_ * layout_->deriv(Ez, {ijk...}, DirectionTag<Direction::X>{});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void BzEq_(Field const& Bz, VecField const& E, Field& Bznew, Indexes const&... ijk) const
    {
        auto const& [Ex, Ey, _] = E();

        if constexpr (dimension == 1)
            Bznew(ijk...)
                = Bz(ijk...) - dt_ * layout_->deriv(Ey, {ijk...}, DirectionTag<Direction::X>{});

        else
            Bznew(ijk...) = Bz(ijk...)
                            - dt_ * layout_->deriv(Ey, {ijk...}, DirectionTag<Direction::X>{})
                            + dt_ * layout_->deriv(Ex, {ijk...}, DirectionTag<Direction::Y>{});
    }
};

} // namespace PHARE::core


#endif
