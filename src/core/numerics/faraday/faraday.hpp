#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP


#include "core/def.hpp"
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
    void operator()(VecField const& B, VecField const& E, VecField& Bnew, double dt) _PHARE_ALL_FN_
    {
        dt_ = dt;
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        layout_.evalOnBox(
            Bxnew, [] _PHARE_ALL_FN_(auto&&... args) { BxEq_(args...); }, Bx, E, Bxnew, *this);
        layout_.evalOnBox(
            Bynew, [] _PHARE_ALL_FN_(auto&&... args) { ByEq_(args...); }, By, E, Bynew, *this);
        layout_.evalOnBox(
            Bznew, [] _PHARE_ALL_FN_(auto&&... args) { BzEq_(args...); }, Bz, E, Bznew, *this);
    }

private:
    GridLayout layout_;
    double dt_;



    template<typename IJK, typename... Args>
    static void BxEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Bx, E, Bxnew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [_, Ey, Ez]          = E();

        if constexpr (dimension == 1)
            Bxnew(ijk) = Bx(ijk);

        if constexpr (dimension == 2)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, ijk);

        if constexpr (dimension == 3)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, ijk)
                         + dt * layout.template deriv<Direction::Z>(Ey, ijk);
    }

    template<typename IJK, typename... Args>
    static void ByEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [By, E, Bynew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [Ex, _, Ez]          = E();

        if constexpr (dimension == 1 || dimension == 2)
            Bynew(ijk) = By(ijk) + dt * layout.template deriv<Direction::X>(Ez, ijk);

        if constexpr (dimension == 3)
            Bynew(ijk) = By(ijk) - dt * layout.template deriv<Direction::Z>(Ex, ijk)
                         + dt * layout.template deriv<Direction::X>(Ez, ijk);
    }

    template<typename IJK, typename... Args>
    static void BzEq_(IJK const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Bz, E, Bznew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [Ex, Ey, _]          = E();

        if constexpr (dimension == 1)
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, ijk);

        else
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, ijk)
                         + dt * layout.template deriv<Direction::Y>(Ex, ijk);
    }
};

} // namespace PHARE::core


#endif
