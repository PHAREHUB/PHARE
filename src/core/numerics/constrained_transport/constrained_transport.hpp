#ifndef PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_HPP
#define PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/index/index.hpp"
#include <iomanip>
#include <tuple>

namespace PHARE::core
{
template<typename GridLayout>
class ConstrainedTransport_ref;

template<typename GridLayout>
class ConstrainedTransport : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    template<typename VecField, typename Fluxes>
    void operator()(VecField& E, Fluxes const& fluxes) const
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - ConstrainedTransport - GridLayout not set, cannot proceed to computation");

        ConstrainedTransport_ref{*this->layout_}(E, fluxes);
    }
};

template<typename GridLayout>
class ConstrainedTransport_ref
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    ConstrainedTransport_ref(GridLayout const& layout)
        : layout_{layout}
    {
    }

    template<typename VecField, typename Fluxes>
    void operator()(VecField& E, Fluxes const& fluxes) const
    {
        auto& Ex = E(Component::X);
        auto& Ey = E(Component::Y);
        auto& Ez = E(Component::Z);


        auto const& By_x = fluxes.B_fx(Component::Y);
        auto const& Bz_x = fluxes.B_fx(Component::Z);

        if constexpr (dimension == 1)
        {
            layout_.evalOnBox(Ey, [&](auto&... args) mutable { EyEq_(Ey, {args...}, Bz_x); });

            layout_.evalOnBox(Ez, [&](auto&... args) mutable { EzEq_(Ez, {args...}, By_x); });
        }
        else if constexpr (dimension >= 2)
        {
            auto const& Bx_y = fluxes.B_fy(Component::X);
            auto const& Bz_y = fluxes.B_fy(Component::Z);

            if constexpr (dimension == 2)
            {
                layout_.evalOnBox(Ex, [&](auto&... args) mutable { ExEq_(Ex, {args...}, Bz_y); });

                layout_.evalOnBox(Ey, [&](auto&... args) mutable { EyEq_(Ey, {args...}, Bz_x); });

                layout_.evalOnBox(Ez,
                                  [&](auto&... args) mutable { EzEq_(Ez, {args...}, By_x, Bx_y); });
            }
            else if constexpr (dimension == 3)
            {
                auto const& Bx_z = fluxes.B_fz(Component::X);
                auto const& By_z = fluxes.B_fz(Component::Y);

                layout_.evalOnBox(Ex,
                                  [&](auto&... args) mutable { ExEq_(Ex, {args...}, Bz_y, By_z); });

                layout_.evalOnBox(Ey,
                                  [&](auto&... args) mutable { EyEq_(Ey, {args...}, Bz_x, Bx_z); });

                layout_.evalOnBox(Ez,
                                  [&](auto&... args) mutable { EzEq_(Ez, {args...}, By_x, Bx_y); });
            }
        }
    }

private:
    GridLayout layout_;

    template<typename Field, typename... Fluxes>
    void ExEq_(Field& Ex, MeshIndex<Field::dimension> index, Fluxes const&... fluxes) const
    {
        auto&& flux_tuple = std::forward_as_tuple(fluxes...);

        if constexpr (dimension >= 2)
        {
            auto& Bz_y = std::get<0>(flux_tuple);

            if constexpr (dimension == 2)
            {
                Ex(index) = -Bz_y(index);
            }
            else if constexpr (dimension == 3)
            {
                auto& By_z = std::get<1>(flux_tuple);

                Ex(index) = 0.25
                            * (-Bz_y(index) - Bz_y(index[0], index[1], index[2] - 1) + By_z(index)
                               + By_z(index[0], index[1] - 1, index[2]));
            }
        }
    }

    template<typename Field, typename... Fluxes>
    void EyEq_(Field& Ey, MeshIndex<Field::dimension> index, Fluxes const&... fluxes) const
    {
        auto&& flux_tuple = std::forward_as_tuple(fluxes...);

        auto& Bz_x = std::get<0>(flux_tuple);

        if constexpr (dimension <= 2)
        {
            Ey(index) = Bz_x(index);
        }
        else if constexpr (dimension == 3)
        {
            auto& Bx_z = std::get<1>(flux_tuple);

            Ey(index) = 0.25
                        * (Bz_x(index) + Bz_x(index[0], index[1], index[2] - 1) - Bx_z(index)
                           - Bx_z(index[0] - 1, index[1], index[2]));
        }
    }

    template<typename Field, typename... Fluxes>
    void EzEq_(Field& Ez, MeshIndex<Field::dimension> index, Fluxes const&... fluxes) const
    {
        auto&& flux_tuple = std::forward_as_tuple(fluxes...);

        auto& By_x = std::get<0>(flux_tuple);

        if constexpr (dimension == 1)
        {
            Ez(index) = -By_x(index);
        }
        else if constexpr (dimension >= 2)
        {
            auto& Bx_y = std::get<1>(flux_tuple);

            if constexpr (dimension == 2)
            {
                Ez(index) = 0.25
                            * (-By_x(index) - By_x(index[0], index[1] - 1) + Bx_y(index)
                               + Bx_y(index[0] - 1, index[1]));
            }
            else if constexpr (dimension == 3)
            {
                Ez(index) = 0.25
                            * (-By_x(index) - By_x(index[0], index[1] - 1, index[2]) + Bx_y(index)
                               + Bx_y(index[0] - 1, index[1], index[2]));
            }
        }
    }
};

} // namespace PHARE::core

#endif
