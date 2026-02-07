#ifndef PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_HPP
#define PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/index/index.hpp"
#include "core/numerics/ohm/ohm.hpp"

#include <iomanip>
#include <strings.h>
#include <tuple>

namespace PHARE::core
{
template<typename GridLayout, bool Resistivity, bool HyperResistivity>
class ConstrainedTransport_ref;

template<typename GridLayout, bool Resistivity, bool HyperResistivity>
class ConstrainedTransport : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    ConstrainedTransport(PHARE::initializer::PHAREDict const& dict)
        : eta_{dict["resistivity"].template to<double>()}
        , nu_{dict["hyper_resistivity"].template to<double>()}
        , hyper_mode_{cppdict::get_value(dict, "hyper_mode", std::string{"constant"}) == "constant"
                          ? HyperMode::constant
                          : HyperMode::spatial}
    {
    }
    template<typename Field, typename VecField, typename Fluxes>
    void operator()(VecField& E, Fluxes const& fluxes, VecField const& J, VecField const& B,
                    Field const& rho) const
    {
        if (!this->hasLayout())
            throw std::runtime_error(
                "Error - ConstrainedTransport - GridLayout not set, cannot proceed to computation");

        ConstrainedTransport_ref<GridLayout, Resistivity, HyperResistivity>{
            *this->layout_, eta_, nu_, hyper_mode_}(E, fluxes, J, B, rho);
    }

private:
    double const eta_;
    double const nu_;
    HyperMode const hyper_mode_;
};

template<typename GridLayout, bool Resistivity, bool HyperResistivity>
class ConstrainedTransport_ref
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    ConstrainedTransport_ref(GridLayout const& layout, double const eta, double const nu,
                             HyperMode const& hyper_mode)
        : layout_{layout}
        , eta_{eta}
        , nu_{nu}
        , hyper_mode_{hyper_mode}
    {
    }

    template<typename Field, typename VecField, typename Fluxes>
    void operator()(VecField& E, Fluxes const& fluxes, VecField const& J, VecField const& B,
                    Field const& rho) const
    {
        auto& Ex = E(Component::X);
        auto& Ey = E(Component::Y);
        auto& Ez = E(Component::Z);

        auto const& Jx = J(Component::X);
        auto const& Jy = J(Component::Y);
        auto const& Jz = J(Component::Z);

        auto const& By_x = fluxes.B_fx(Component::Y);
        auto const& Bz_x = fluxes.B_fx(Component::Z);

        if constexpr (dimension == 1)
        {
            layout_.evalOnBox(Ey, [&](auto&... args) mutable { EyEq_(Ey, {args...}, Bz_x); });

            layout_.evalOnBox(Ez, [&](auto&... args) mutable { EzEq_(Ez, {args...}, By_x); });

            if constexpr (Resistivity)
            {
                layout_.evalOnBox(
                    Ey, [&](auto&... args) mutable { resistive_contribution_(Ey, Jy, {args...}); });
                layout_.evalOnBox(
                    Ez, [&](auto&... args) mutable { resistive_contribution_(Ez, Jz, {args...}); });
            }

            if constexpr (HyperResistivity)
            {
                layout_.evalOnBox(Ey, [&](auto&... args) mutable {
                    hyperresistive_contribution_<Component::Y>(Ey, Jy, B, rho, {args...});
                });
                layout_.evalOnBox(Ez, [&](auto&... args) mutable {
                    hyperresistive_contribution_<Component::Z>(Ez, Jz, B, rho, {args...});
                });
            }
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

            if constexpr (Resistivity)
            {
                layout_.evalOnBox(
                    Ex, [&](auto&... args) mutable { resistive_contribution_(Ex, Jx, {args...}); });
                layout_.evalOnBox(
                    Ey, [&](auto&... args) mutable { resistive_contribution_(Ey, Jy, {args...}); });
                layout_.evalOnBox(
                    Ez, [&](auto&... args) mutable { resistive_contribution_(Ez, Jz, {args...}); });
            }

            if constexpr (HyperResistivity)
            {
                layout_.evalOnBox(Ex, [&](auto&... args) mutable {
                    hyperresistive_contribution_<Component::X>(Ex, Jx, B, rho, {args...});
                });
                layout_.evalOnBox(Ey, [&](auto&... args) mutable {
                    hyperresistive_contribution_<Component::Y>(Ey, Jy, B, rho, {args...});
                });
                layout_.evalOnBox(Ez, [&](auto&... args) mutable {
                    hyperresistive_contribution_<Component::Z>(Ez, Jz, B, rho, {args...});
                });
            }
        }
    }

private:
    GridLayout layout_;
    double const eta_;
    double const nu_;
    HyperMode const hyper_mode_;

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

    template<typename Field>
    void resistive_contribution_(Field& E, Field const& J, MeshIndex<Field::dimension> index) const
    {
        E(index) += eta_ * J(index);
    }

    template<auto component, typename Field, typename VecField>
    void hyperresistive_contribution_(Field& E, Field const& J, VecField const& B, Field const& rho,
                                      MeshIndex<Field::dimension> index) const
    {
        if (hyper_mode_ == HyperMode::constant)
            return constant_hyperresistive_<component>(E, J, index);
        else if (hyper_mode_ == HyperMode::spatial)
            return spatial_hyperresistive_<component>(E, J, B, rho, index);
        else
            throw std::runtime_error("Error - Ohm - unknown hyper_mode");
    }

    template<auto component, typename Field>
    void constant_hyperresistive_(Field& E, Field const& J, MeshIndex<Field::dimension> index) const
    {
        E(index) -= nu_ * layout_.laplacian(J, index);
    }

    template<auto component, typename Field, typename VecField>
    void spatial_hyperresistive_(Field& E, Field const& J, VecField const& B, Field const& rho,
                                 MeshIndex<Field::dimension> index) const
    {
        auto minMeshSize = [&]() {
            auto const meshSize = layout_.meshSize();
            if constexpr (Field::dimension == 1)
                return meshSize[0];
            else if constexpr (Field::dimension == 2)
                return std::min({meshSize[0], meshSize[1]});
            else
                return std::min({meshSize[0], meshSize[1], meshSize[2]});
        }();

        auto computeHR = [&](auto BxProj, auto ByProj, auto BzProj, auto rhoProj) {
            auto const BxOnE = GridLayout::project(B(Component::X), index, BxProj);
            auto const ByOnE = GridLayout::project(B(Component::Y), index, ByProj);
            auto const BzOnE = GridLayout::project(B(Component::Z), index, BzProj);
            auto const nOnE  = GridLayout::project(rho, index, rhoProj);
            auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);
            E(index)
                -= nu_ * layout_.laplacian(J, index) * minMeshSize * minMeshSize * (b / nOnE + 1);
        };

        if constexpr (component == Component::X)
        {
            return computeHR(GridLayout::BxToEx(), GridLayout::ByToEx(), GridLayout::BzToEx(),
                             GridLayout::cellCenterToEdgeX());
        }
        if constexpr (component == Component::Y)
        {
            return computeHR(GridLayout::BxToEy(), GridLayout::ByToEy(), GridLayout::BzToEy(),
                             GridLayout::cellCenterToEdgeY());
        }
        if constexpr (component == Component::Z)
        {
            return computeHR(GridLayout::BxToEz(), GridLayout::ByToEz(), GridLayout::BzToEz(),
                             GridLayout::cellCenterToEdgeZ());
        }
    }
};

} // namespace PHARE::core

#endif
