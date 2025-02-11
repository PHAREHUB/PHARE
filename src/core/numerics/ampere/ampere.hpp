#ifndef PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP
#define PHARE_CORE_NUMERICS_AMPERE_AMPERE_HPP

#include <cstddef>
#include <iostream>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"


namespace PHARE::core
{
template<typename GridLayout>
class Ampere_ref;

template<typename GridLayout>
class Ampere : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField>
    void operator()(VecField const& B, VecField& J)
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - Ampere - GridLayout not set, cannot proceed to calculate ampere()");

        Ampere_ref{*this->layout_}(B, J);
    }
};

template<typename GridLayout>
class Ampere_ref
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    Ampere_ref(GridLayout const& layout)
        : layout_{layout}
    {
    }

    template<typename VecField>
    void operator()(VecField const& B, VecField& J) const
    {
        // can't use structured bindings because
        //   "reference to local binding declared in enclosing function"
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        layout_.evalOnBox(Jx, [&](auto&... args) mutable { JxEq_(Jx, B, args...); });
        layout_.evalOnBox(Jy, [&](auto&... args) mutable { JyEq_(Jy, B, args...); });
        layout_.evalOnBox(Jz, [&](auto&... args) mutable { JzEq_(Jz, B, args...); });
    }


private:
    GridLayout layout_;

    template<typename VecField, typename Field, typename... Indexes>
    void JxEq_(Field& Jx, VecField const& B, Indexes const&... ijk) const
    {
        auto const& [_, By, Bz] = B();

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
