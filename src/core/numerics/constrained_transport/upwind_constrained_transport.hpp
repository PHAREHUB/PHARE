#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_HPP


#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/numerics/constrained_transport/upwind_edge_b.hpp"
#include "core/utilities/index/index.hpp"

#include <cmath>

namespace PHARE::core
{
using UpwindConstrainedTransportInfo = OhmInfo;

template<typename GridLayout, template<typename> typename Reconstruction, bool Hall>
class UpwindConstrainedTransport : UpwindConstrainedTransportInfo
{
    using Super                     = UpwindConstrainedTransportInfo;
    using Reconstruction_t          = Reconstruction<GridLayout>;
    using EdgeB_t                   = UpwindEdgeB<GridLayout, Reconstruction>;
    constexpr static auto dimension = GridLayout::dimension;
    using Super::hyper_mode;

public:
    using Info_t = Super;

    UpwindConstrainedTransport(UpwindConstrainedTransportInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
        , is_dissipative_{info.isResistive() || info.isHyperResistive()}
    {
    }

    void operator()(auto& ct_state, auto const& dissipative_electric_state, auto& mhd_state) const
    {
        auto& E       = mhd_state.E;
        auto const& B = mhd_state.B;

        auto& Ex = E(Component::X);
        auto& Ey = E(Component::Y);
        auto& Ez = E(Component::Z);

        layout_.evalOnBox(Ex, [&](auto&... args) { ExEq_(ct_state, Ex, B, {args...}); });
        layout_.evalOnBox(Ey, [&](auto&... args) { EyEq_(ct_state, Ey, B, {args...}); });
        layout_.evalOnBox(Ez, [&](auto&... args) { EzEq_(ct_state, Ez, B, {args...}); });

        if (is_dissipative_)
        {
            auto const& Ediss = dissipative_electric_state.E();
            for (auto c : {Component::X, Component::Y, Component::Z})
                layout_.evalOnBox(E(c), [&](auto&... args) { E(c)(args...) += Ediss(c)(args...); });
        }
    }

private:
    void ExEq_(auto& ct_state, auto& Ex, auto const& B, MeshIndex<dimension> idx) const
    {
        if constexpr (dimension == 2)
        {
            auto const [aL, aR, dL, dR, BzL, BzR]
                = EdgeB_t::template states<Component::X, Component::Z>(ct_state, B, idx);

            auto FL = BzL * ct_state.vt_y(Component::Y)(idx)
                      - B(Component::Y)(idx) * ct_state.vt_y(Component::Z)(idx);
            auto FR = BzR * ct_state.vt_y(Component::Y)(idx)
                      - B(Component::Y)(idx) * ct_state.vt_y(Component::Z)(idx);

            Ex(idx) = -(aL * FL + aR * FR) + (dR * BzR - dL * BzL);

            if constexpr (Hall)
            {
                auto invRho  = 1.0 / ct_state.template getRhot<Direction::Y>()(idx);
                auto JxB_x_L = ct_state.template getJt<Direction::Y>()(Component::Y)(idx)*BzL
                               - ct_state.template getJt<Direction::Y>()(Component::Z)(idx)*B(
                                   Component::Y)(idx);
                auto JxB_x_R = ct_state.template getJt<Direction::Y>()(Component::Y)(idx)*BzR
                               - ct_state.template getJt<Direction::Y>()(Component::Z)(idx)*B(
                                   Component::Y)(idx);

                auto HallL = -JxB_x_L * invRho;
                auto HallR = -JxB_x_R * invRho;

                auto F_Bz_y = aL * HallL + aR * HallR;

                Ex(idx) += -F_Bz_y;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto const [aS, aN, dS, dN, BzS, BzN]
                = EdgeB_t::template states<Component::X, Component::Z>(ct_state, B, idx);
            auto const [aB, aT, dB, dT, ByB, ByT]
                = EdgeB_t::template states<Component::X, Component::Y>(ct_state, B, idx);

            auto [vyS, vyN] = Reconstruction_t::template reconstruct<Direction::Y>(
                ct_state.vt_z(Component::Y), idx);
            auto [vzB, vzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                ct_state.vt_y(Component::Z), idx);

            Ex(idx) = (aB * vzB * ByB + aT * vzT * ByT) - (aS * vyS * BzS + aN * vyN * BzN)
                      - (dT * ByT - dB * ByB) + (dN * BzN - dS * BzS);

            if constexpr (Hall)
            {
                auto [jyS, jyN] = Reconstruction_t::template reconstruct<Direction::Y>(
                    ct_state.template getJt<Direction::Z>()(Component::Y), idx);
                auto [jzB, jzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getJt<Direction::Y>()(Component::Z), idx);

                auto [rhoS, rhoN] = Reconstruction_t::template reconstruct<Direction::Y>(
                    ct_state.template getRhot<Direction::Z>(), idx);
                auto [rhoB, rhoT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getRhot<Direction::Y>(), idx);

                Ex(idx) += -(aB * jzB * ByB / rhoB + aT * jzT * ByT / rhoT)
                           + (aS * jyS * BzS / rhoS + aN * jyN * BzN / rhoN);
            }
        }
    }

    void EyEq_(auto& ct_state, auto& Ey, auto const& B, MeshIndex<dimension> idx) const
    {
        if constexpr (dimension <= 2)
        {
            auto const [aL, aR, dL, dR, BzL, BzR]
                = EdgeB_t::template states<Component::Y, Component::Z>(ct_state, B, idx);

            auto FL = BzL * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Z)(idx);
            auto FR = BzR * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Z)(idx);

            Ey(idx) = (aL * FL + aR * FR) - (dR * BzR - dL * BzL);

            if constexpr (Hall)
            {
                auto invRho  = 1.0 / ct_state.template getRhot<Direction::X>()(idx);
                auto JxB_y_L = ct_state.template getJt<Direction::X>()(Component::Z)(idx)*B(
                                   Component::X)(idx)
                               - ct_state.template getJt<Direction::X>()(Component::X)(idx)*BzL;
                auto JxB_y_R = ct_state.template getJt<Direction::X>()(Component::Z)(idx)*B(
                                   Component::X)(idx)
                               - ct_state.template getJt<Direction::X>()(Component::X)(idx)*BzR;

                auto HallL = JxB_y_L * invRho;
                auto HallR = JxB_y_R * invRho;

                auto F_Bz_x = aL * HallL + aR * HallR;

                Ey(idx) += F_Bz_x;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto const [aW, aE, dW, dE, BzW, BzE]
                = EdgeB_t::template states<Component::Y, Component::Z>(ct_state, B, idx);
            auto const [aB, aT, dB, dT, BxB, BxT]
                = EdgeB_t::template states<Component::Y, Component::X>(ct_state, B, idx);

            auto [vxW, vxE] = Reconstruction_t::template reconstruct<Direction::X>(
                ct_state.vt_z(Component::X), idx);
            auto [vzB, vzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                ct_state.vt_x(Component::Z), idx);

            Ey(idx) = (aW * vxW * BzW + aE * vxE * BzE) - (aB * vzB * BxB + aT * vzT * BxT)
                      - (dE * BzE - dW * BzW) + (dT * BxT - dB * BxB);

            if constexpr (Hall)
            {
                auto [jxW, jxE] = Reconstruction_t::template reconstruct<Direction::X>(
                    ct_state.template getJt<Direction::Z>()(Component::X), idx);
                auto [jzB, jzT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getJt<Direction::X>()(Component::Z), idx);
                auto [rhoW, rhoE] = Reconstruction_t::template reconstruct<Direction::X>(
                    ct_state.template getRhot<Direction::Z>(), idx);
                auto [rhoB, rhoT] = Reconstruction_t::template reconstruct<Direction::Z>(
                    ct_state.template getRhot<Direction::X>(), idx);
                Ey(idx) += -(aW * jxW * BzW / rhoW + aE * jxE * BzE / rhoE)
                           + (aB * jzB * BxB / rhoB + aT * jzT * BxT / rhoT);
            }
        }
    }

    void EzEq_(auto& ct_state, auto& Ez, auto const& B, MeshIndex<dimension> idx) const
    {
        if constexpr (dimension == 1)
        {
            auto const [aL, aR, dL, dR, ByL, ByR]
                = EdgeB_t::template states<Component::Z, Component::Y>(ct_state, B, idx);

            auto FL = ByL * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Y)(idx);
            auto FR = ByR * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Y)(idx);

            Ez(idx) = -(aL * FL + aR * FR) + (dR * ByR - dL * ByL);

            if constexpr (Hall)
            {
                auto invRho  = 1.0 / ct_state.template getRhot<Direction::X>()(idx);
                auto JxB_z_L = ct_state.template getJt<Direction::X>()(Component::X)(idx)*ByL
                               - ct_state.template getJt<Direction::X>()(Component::Y)(idx)*B(
                                   Component::X)(idx);
                auto JxB_z_R = ct_state.template getJt<Direction::X>()(Component::X)(idx)*ByR
                               - ct_state.template getJt<Direction::X>()(Component::Y)(idx)*B(
                                   Component::X)(idx);

                auto HallL = -JxB_z_L * invRho;
                auto HallR = -JxB_z_R * invRho;

                auto F_By_x = aL * HallL + aR * HallR;

                Ez(idx) += -F_By_x;
            }
        }
        else if constexpr (dimension >= 2)
        {
            auto const [aW, aE, dW, dE, ByW, ByE]
                = EdgeB_t::template states<Component::Z, Component::Y>(ct_state, B, idx);
            auto const [aS, aN, dS, dN, BxS, BxN]
                = EdgeB_t::template states<Component::Z, Component::X>(ct_state, B, idx);

            auto [vyS, vyN] = Reconstruction_t::template reconstruct<Direction::Y>(
                ct_state.vt_x(Component::Y), idx);
            auto [vxW, vxE] = Reconstruction_t::template reconstruct<Direction::X>(
                ct_state.vt_y(Component::X), idx);

            Ez(idx) = -(aW * vxW * ByW + aE * vxE * ByE) + (aS * vyS * BxS + aN * vyN * BxN)
                      + (dE * ByE - dW * ByW) - (dN * BxN - dS * BxS);

            if constexpr (Hall)
            {
                auto [jyS, jyN] = Reconstruction_t::template reconstruct<Direction::Y>(
                    ct_state.template getJt<Direction::X>()(Component::Y), idx);
                auto [jxW, jxE] = Reconstruction_t::template reconstruct<Direction::X>(
                    ct_state.template getJt<Direction::Y>()(Component::X), idx);

                auto [rhoS, rhoN] = Reconstruction_t::template reconstruct<Direction::Y>(
                    ct_state.template getRhot<Direction::X>(), idx);
                auto [rhoW, rhoE] = Reconstruction_t::template reconstruct<Direction::X>(
                    ct_state.template getRhot<Direction::Y>(), idx);

                Ez(idx) += (aW * jxW * ByW / rhoW + aE * jxE * ByE / rhoE)
                           - (aS * jyS * BxS / rhoS + aN * jyN * BxN / rhoN);
            }
        }
    }

    GridLayout layout_;
    bool const is_dissipative_;
};
} // namespace PHARE::core

#endif
