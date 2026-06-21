#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"


namespace PHARE::core
{



template<typename GridLayout>
class Faraday
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Faraday(GridLayout const& layout)
        : layout_{layout}
    {
    }


    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt)
    {
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

        layout_.evalOnBox(
            Bxnew, [](auto&&... args) { BxEq_(args...); }, Bx, E, Bxnew, layout_, dt_);
        layout_.evalOnBox(
            Bynew, [](auto&&... args) { ByEq_(args...); }, By, E, Bynew, layout_, dt_);
        layout_.evalOnBox(
            Bznew, [](auto&&... args) { BzEq_(args...); }, Bz, E, Bznew, layout_, dt_);
    }

private:
    double dt_;
    GridLayout layout_;

    static void BxEq_(auto const& ijk, auto&&... args)
    {
        auto const& [Bx, E, Bxnew, layout, dt] = std::forward_as_tuple(args...);
        auto const& [_, Ey, Ez]                = E();

        if constexpr (dimension == 1)
            Bxnew(ijk) = Bx(ijk);

        if constexpr (dimension == 2)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, ijk);

        if constexpr (dimension == 3)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, ijk)
                         + dt * layout.template deriv<Direction::Z>(Ey, ijk);
    }


    static void ByEq_(auto const& ijk, auto&&... args)
    {
        auto const& [By, E, Bynew, layout, dt] = std::forward_as_tuple(args...);
        auto const& [Ex, _, Ez]                = E();

        if constexpr (dimension == 1 || dimension == 2)
            Bynew(ijk) = By(ijk) + dt * layout.template deriv<Direction::X>(Ez, ijk);

        if constexpr (dimension == 3)
            Bynew(ijk) = By(ijk) - dt * layout.template deriv<Direction::Z>(Ex, ijk)
                         + dt * layout.template deriv<Direction::X>(Ez, ijk);
    }


    static void BzEq_(auto const& ijk, auto&&... args)
    {
        auto const& [Bz, E, Bznew, layout, dt] = std::forward_as_tuple(args...);
        auto const& [Ex, Ey, _]                = E();

        if constexpr (dimension == 1)
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, ijk);

        else
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, ijk)
                         + dt * layout.template deriv<Direction::Y>(Ex, ijk);
    }
};


} // namespace PHARE::core


#endif
