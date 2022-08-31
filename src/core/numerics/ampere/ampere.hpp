#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include <cstddef>
#include <iostream>

#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"


namespace PHARE::core
{
template<typename GridLayout>
class Ampere : public LayoutHolder<GridLayout>
{
    auto static constexpr dimension = GridLayout::dimension;

public:
    template<typename VecField>
    void operator()(VecField const& B_, VecField& J_)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ampere - GridLayout not set, cannot proceed to calculate ampere()");

        auto B = B_.view();
        auto J = J_.view();

        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx     = J(Component::X);
        auto& Jy     = J(Component::Y);
        auto& Jz     = J(Component::Z);
        auto& layout = *this->layout_;

        layout.evalOnBox(
            Jx, [] _PHARE_ALL_FN_(auto&... args) { JxEq_(args...); }, Jx, B, layout);
        layout.evalOnBox(
            Jy, [] _PHARE_ALL_FN_(auto&... args) { JyEq_(args...); }, Jy, B, layout);
        layout.evalOnBox(
            Jz, [] _PHARE_ALL_FN_(auto&... args) { JzEq_(args...); }, Jz, B, layout);
    }

private:
    template<typename... IJK, typename... Args>
    static void JxEq_(std::tuple<IJK...> const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Jx, B, layout] = std::forward_as_tuple(args...);
        auto const& [_, By, Bz]     = B();

        if constexpr (dimension == 2)
            Jx(ijk) = layout.template deriv<Direction::Y>(Bz, {ijk});

        if constexpr (dimension == 3)
            Jx(ijk) = layout.template deriv<Direction::Y>(Bz, {ijk})
                      - layout.template deriv<Direction::Z>(By, {ijk});
    }

    template<typename... IJK, typename... Args>
    static void JyEq_(std::tuple<IJK...> const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Jy, B, layout] = std::forward_as_tuple(args...);
        auto const& [Bx, By, Bz]    = B();

        if constexpr (dimension == 1 || dimension == 2)
            Jy(ijk) = -layout.template deriv<Direction::X>(Bz, {ijk});

        if constexpr (dimension == 3)
            Jy(ijk) = layout.template deriv<Direction::Z>(Bx, {ijk})
                      - layout.template deriv<Direction::X>(Bz, {ijk});
    }

    template<typename... IJK, typename... Args>
    static void JzEq_(std::tuple<IJK...> const& ijk, Args&&... args) _PHARE_ALL_FN_
    {
        auto const& [Jz, B, layout] = std::forward_as_tuple(args...);
        auto const& [Bx, By, Bz]    = B();

        if constexpr (dimension == 1)
            Jz(ijk) = layout.template deriv<Direction::X>(By, {ijk});

        else
            Jz(ijk) = layout.template deriv<Direction::X>(By, {ijk})
                      - layout.template deriv<Direction::Y>(Bx, {ijk});
    }
};

} // namespace PHARE::core
#endif
