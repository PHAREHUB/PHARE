#ifndef PHARE_FARADAY_HPP
#define PHARE_FARADAY_HPP

#include <cstddef>

#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"


namespace PHARE::core
{
template<typename GridLayout>
class Faraday_ref;

template<typename GridLayout>
class Faraday : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    template<typename VecField>
    void operator()(VecField const& B_, VecField const& E_, VecField& Bnew_,
                    double dt) _PHARE_ALL_FN_
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Faraday - GridLayout not set, cannot proceed to calculate faraday()");

        if (!(B_.isUsable() && E_.isUsable() && Bnew_.isUsable()))
            throw std::runtime_error("Error - Faraday - not all VecField parameters are usable");

        auto E    = E_.view();
        auto B    = B_.view();
        auto Bnew = Bnew_.view();

        Faraday_ref{*this->layout_, dt}(B, E, Bnew);
    }
};

template<typename GridLayout>
class Faraday_ref
{
    constexpr static auto dimension = GridLayout::dimension;
    using This                      = Faraday_ref<GridLayout>;

public:
    Faraday_ref(GridLayout& layout_, double dt_)
        : layout{layout_}
        , dt{dt_}
    {
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField const& E, VecField& Bnew) _PHARE_ALL_FN_
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto& Bxnew = Bnew(Component::X);
        auto& Bynew = Bnew(Component::Y);
        auto& Bznew = Bnew(Component::Z);

        layout.evalOnBox(
            Bxnew, [] _PHARE_ALL_FN_(auto&... args) { BxEq_(args...); }, Bx, E, Bxnew, *this);
        layout.evalOnBox(
            Bynew, [] _PHARE_ALL_FN_(auto&... args) { ByEq_(args...); }, By, E, Bynew, *this);
        layout.evalOnBox(
            Bznew, [] _PHARE_ALL_FN_(auto&... args) { BzEq_(args...); }, Bz, E, Bznew, *this);
    }

private:
    GridLayout layout;
    double dt;

    template<typename... IJK, typename... Args>
    static void BxEq_(std::tuple<IJK...> const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Bx, E, Bxnew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [_, Ey, Ez]          = E();

        if constexpr (dimension == 1)
            Bxnew(ijk) = Bx(ijk);

        if constexpr (dimension == 2)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, {ijk});

        if constexpr (dimension == 3)
            Bxnew(ijk) = Bx(ijk) - dt * layout.template deriv<Direction::Y>(Ez, {ijk})
                         + dt * layout.template deriv<Direction::Z>(Ey, {ijk});
    }

    template<typename... IJK, typename... Args>
    static void ByEq_(std::tuple<IJK...> const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [By, E, Bynew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [Ex, _, Ez]          = E();

        if constexpr (dimension == 1 || dimension == 2)
            Bynew(ijk) = By(ijk) + dt * layout.template deriv<Direction::X>(Ez, {ijk});

        if constexpr (dimension == 3)
            Bynew(ijk) = By(ijk) - dt * layout.template deriv<Direction::Z>(Ex, {ijk})
                         + dt * layout.template deriv<Direction::X>(Ez, {ijk});
    }

    template<typename... IJK, typename... Args>
    static void BzEq_(std::tuple<IJK...> const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Bz, E, Bznew, self] = std::forward_as_tuple(args...);
        auto const& [layout, dt]         = self;
        auto const& [Ex, Ey, _]          = E();

        if constexpr (dimension == 1)
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, {ijk});

        else
            Bznew(ijk) = Bz(ijk) - dt * layout.template deriv<Direction::X>(Ey, {ijk})
                         + dt * layout.template deriv<Direction::Y>(Ex, {ijk});
    }
};

} // namespace PHARE::core


#endif
