#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"


namespace PHARE::core
{
template<typename GridLayout>
class Faraday : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
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

        layout_->evalOnBox(Bxnew, [&](auto&... args) mutable { BxEq_(Bx, E, Bxnew, args...); });
        layout_->evalOnBox(Bynew, [&](auto&... args) mutable { ByEq_(By, E, Bynew, args...); });
        layout_->evalOnBox(Bznew, [&](auto&... args) mutable { BzEq_(Bz, E, Bznew, args...); });
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
            Bxnew(ijk...) = Bx(ijk...) - dt_ * layout_->template deriv<Direction::Y>(Ez, {ijk...});

        if constexpr (dimension == 3)
            Bxnew(ijk...) = Bx(ijk...) - dt_ * layout_->template deriv<Direction::Y>(Ez, {ijk...})
                            + dt_ * layout_->template deriv<Direction::Z>(Ey, {ijk...});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void ByEq_(Field const& By, VecField const& E, Field& Bynew, Indexes const&... ijk) const
    {
        auto const& [Ex, _, Ez] = E();

        if constexpr (dimension == 1 || dimension == 2)
            Bynew(ijk...) = By(ijk...) + dt_ * layout_->template deriv<Direction::X>(Ez, {ijk...});

        if constexpr (dimension == 3)
            Bynew(ijk...) = By(ijk...) - dt_ * layout_->template deriv<Direction::Z>(Ex, {ijk...})
                            + dt_ * layout_->template deriv<Direction::X>(Ez, {ijk...});
    }

    template<typename VecField, typename Field, typename... Indexes>
    void BzEq_(Field const& Bz, VecField const& E, Field& Bznew, Indexes const&... ijk) const
    {
        auto const& [Ex, Ey, _] = E();

        if constexpr (dimension == 1)
            Bznew(ijk...) = Bz(ijk...) - dt_ * layout_->template deriv<Direction::X>(Ey, {ijk...});

        else
            Bznew(ijk...) = Bz(ijk...) - dt_ * layout_->template deriv<Direction::X>(Ey, {ijk...})
                            + dt_ * layout_->template deriv<Direction::Y>(Ex, {ijk...});
    }
};

} // namespace PHARE::core


#endif
