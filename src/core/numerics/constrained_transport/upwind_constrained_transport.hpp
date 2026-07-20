#ifndef PHARE_UPWIND_CONSTRAINED_TRANSPORT_HPP
#define PHARE_UPWIND_CONSTRAINED_TRANSPORT_HPP

#include "core/numerics/reconstructions/constant.hpp"
#include "core/def.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/constrained_transport/upwind_constrained_transport_utils.hpp"


#include <cmath>

namespace PHARE::core
{
using UpwindConstrainedTransportInfo = OhmInfo;

template<typename GridLayout, template<typename> typename Reconstruction, bool Hall,
         bool Resistivity, bool HyperResistivity>
class UpwindConstrainedTransport : UpwindConstrainedTransportInfo
{
    using Super                     = UpwindConstrainedTransportInfo;
    using Reconstruction_t          = Reconstruction<GridLayout>;
    constexpr static auto dimension = GridLayout::dimension;
    using Super::hyper_mode;

public:
    using Info_t = Super;

    UpwindConstrainedTransport(UpwindConstrainedTransportInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
    {
    }

    void operator()(auto& ct_state, auto& mhd_state) const
    {
        auto& E       = mhd_state.E;
        auto const& B = mhd_state.B;

        auto& Ex = E(Component::X);
        auto& Ey = E(Component::Y);
        auto& Ez = E(Component::Z);

        layout_.evalOnBox(Ex, [&](auto&... args) mutable { ExEq_(ct_state, Ex, B, {args...}); });
        layout_.evalOnBox(Ey, [&](auto&... args) mutable { EyEq_(ct_state, Ey, B, {args...}); });
        layout_.evalOnBox(Ez, [&](auto&... args) mutable { EzEq_(ct_state, Ez, B, {args...}); });

        if constexpr (Resistivity || HyperResistivity)
        {
            auto const& J = mhd_state.J;

            auto& Jx = J(Component::X);
            auto& Jy = J(Component::Y);
            auto& Jz = J(Component::Z);

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
                auto const& rho = mhd_state.rho;

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

    /**
     * @brief Inner-boundary degraded-E correction (second pass).
     *
     * Recompute E with first-order (piecewise-constant) reconstruction and an ideal Ohm's law at
     * the edges near an under-resolved inner boundary, consistently with the degraded hydro
     * fluxes. The edge lists are keyed by the E-component centering and are empty on resolved
     * levels, so this is a no-op there.
     */
    template<typename MeshData>
    void degrade_E_near_inner_boundary(auto& ct_state, auto& mhd_state,
                                       MeshData const& meshData) const
    {
        auto& E       = mhd_state.E;
        auto const& B = mhd_state.B;

        using Constant = ConstantReconstruction<GridLayout>;

        auto degradeComp = [&](auto& Ecomp, auto&& recompute) {
            auto const centering = layout_.centering(Ecomp.physicalQuantity());
            auto const& degraded = meshData.getDegradedElemsFromCentering(centering);
            for (auto const& elem : degraded)
                recompute(elem.index.template as<std::size_t>());
        };

        degradeComp(E(Component::X), [&](auto const& idx) {
            ExEq_<Constant, true>(ct_state, E(Component::X), B, idx);
        });
        degradeComp(E(Component::Y), [&](auto const& idx) {
            EyEq_<Constant, true>(ct_state, E(Component::Y), B, idx);
        });
        degradeComp(E(Component::Z), [&](auto const& idx) {
            EzEq_<Constant, true>(ct_state, E(Component::Z), B, idx);
        });
    }

private:
    template<typename Rec = Reconstruction_t, bool IdealOnly = false>
    void ExEq_(auto& ct_state, auto& Ex, auto const& B, MeshIndex<dimension> idx) const
    {
        if constexpr (dimension == 2)
        {
            auto [BzL, BzR]
                = Rec::template reconstruct<Direction::Y>(B(Component::Z), idx);

            auto FL = BzL * ct_state.vt_y(Component::Y)(idx)
                      - B(Component::Y)(idx) * ct_state.vt_y(Component::Z)(idx);
            auto FR = BzR * ct_state.vt_y(Component::Y)(idx)
                      - B(Component::Y)(idx) * ct_state.vt_y(Component::Z)(idx);

            Ex(idx) = -(ct_state.aL_y(idx) * FL + ct_state.aR_y(idx) * FR)
                      + (ct_state.dR_y(idx) * BzR - ct_state.dL_y(idx) * BzL);

            if constexpr (Hall && !IdealOnly)
            {
                auto invRho  = 1.0 / ct_state.rhot_y(idx);
                auto JxB_x_L = ct_state.jt_y(Component::Y)(idx) * BzL
                               - ct_state.jt_y(Component::Z)(idx) * B(Component::Y)(idx);
                auto JxB_x_R = ct_state.jt_y(Component::Y)(idx) * BzR
                               - ct_state.jt_y(Component::Z)(idx) * B(Component::Y)(idx);

                auto HallL = -JxB_x_L * invRho;
                auto HallR = -JxB_x_R * invRho;

                auto F_Bz_y = ct_state.aL_y(idx) * HallL + ct_state.aR_y(idx) * HallR;

                Ex(idx) += -F_Bz_y;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto aS = 0.5
                      * (ct_state.aL_y(idx)
                         + ct_state.aL_y(layout_.template previous<Direction::Z>(idx)));
            auto aN = 0.5
                      * (ct_state.aR_y(idx)
                         + ct_state.aR_y(layout_.template previous<Direction::Z>(idx)));
            auto aB = 0.5
                      * (ct_state.aL_z(idx)
                         + ct_state.aL_z(layout_.template previous<Direction::Y>(idx)));
            auto aT = 0.5
                      * (ct_state.aR_z(idx)
                         + ct_state.aR_z(layout_.template previous<Direction::Y>(idx)));
            auto dS = 0.5
                      * (ct_state.dL_y(idx)
                         + ct_state.dL_y(layout_.template previous<Direction::Z>(idx)));
            auto dN = 0.5
                      * (ct_state.dR_y(idx)
                         + ct_state.dR_y(layout_.template previous<Direction::Z>(idx)));
            auto dB = 0.5
                      * (ct_state.dL_z(idx)
                         + ct_state.dL_z(layout_.template previous<Direction::Y>(idx)));
            auto dT = 0.5
                      * (ct_state.dR_z(idx)
                         + ct_state.dR_z(layout_.template previous<Direction::Y>(idx)));

            auto [vyS, vyN] = Rec::template reconstruct<Direction::Y>(
                ct_state.vt_z(Component::Y), idx);
            auto [vzB, vzT] = Rec::template reconstruct<Direction::Z>(
                ct_state.vt_y(Component::Z), idx);

            auto [BzS, BzN]
                = Rec::template reconstruct<Direction::Y>(B(Component::Z), idx);
            auto [ByB, ByT]
                = Rec::template reconstruct<Direction::Z>(B(Component::Y), idx);

            Ex(idx) = (aB * vzB * ByB + aT * vzT * ByT) - (aS * vyS * BzS + aN * vyN * BzN)
                      - (dT * ByT - dB * ByB) + (dN * BzN - dS * BzS);

            if constexpr (Hall && !IdealOnly)
            {
                auto [jyS, jyN] = Rec::template reconstruct<Direction::Y>(
                    ct_state.jt_z(Component::Y), idx);
                auto [jzB, jzT] = Rec::template reconstruct<Direction::Z>(
                    ct_state.jt_y(Component::Z), idx);

                auto [rhoS, rhoN]
                    = Rec::template reconstruct<Direction::Y>(ct_state.rhot_z, idx);
                auto [rhoB, rhoT]
                    = Rec::template reconstruct<Direction::Z>(ct_state.rhot_y, idx);

                Ex(idx) += -(aB * jzB * ByB / rhoB + aT * jzT * ByT / rhoT)
                           + (aS * jyS * BzS / rhoS + aN * jyN * BzN / rhoN);
            }
        }
    }

    template<typename Rec = Reconstruction_t, bool IdealOnly = false>
    void EyEq_(auto& ct_state, auto& Ey, auto const& B, MeshIndex<dimension> idx) const
    {
        if constexpr (dimension <= 2)
        {
            auto [BzL, BzR]
                = Rec::template reconstruct<Direction::X>(B(Component::Z), idx);

            auto FL = BzL * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Z)(idx);
            auto FR = BzR * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Z)(idx);

            Ey(idx) = (ct_state.aL_x(idx) * FL + ct_state.aR_x(idx) * FR)
                      - (ct_state.dR_x(idx) * BzR - ct_state.dL_x(idx) * BzL);

            if constexpr (Hall && !IdealOnly)
            {
                auto invRho  = 1.0 / ct_state.rhot_x(idx);
                auto JxB_y_L = ct_state.jt_x(Component::Z)(idx) * B(Component::X)(idx)
                               - ct_state.jt_x(Component::X)(idx) * BzL;
                auto JxB_y_R = ct_state.jt_x(Component::Z)(idx) * B(Component::X)(idx)
                               - ct_state.jt_x(Component::X)(idx) * BzR;

                auto HallL = JxB_y_L * invRho;
                auto HallR = JxB_y_R * invRho;

                auto F_Bz_x = ct_state.aL_x(idx) * HallL + ct_state.aR_x(idx) * HallR;

                Ey(idx) += F_Bz_x;
            }
        }
        else if constexpr (dimension == 3)
        {
            auto aW = 0.5
                      * (ct_state.aL_x(idx)
                         + ct_state.aL_x(layout_.template previous<Direction::Z>(idx)));
            auto aE = 0.5
                      * (ct_state.aR_x(idx)
                         + ct_state.aR_x(layout_.template previous<Direction::Z>(idx)));
            auto aB = 0.5
                      * (ct_state.aL_z(idx)
                         + ct_state.aL_z(layout_.template previous<Direction::X>(idx)));
            auto aT = 0.5
                      * (ct_state.aR_z(idx)
                         + ct_state.aR_z(layout_.template previous<Direction::X>(idx)));
            auto dW = 0.5
                      * (ct_state.dL_x(idx)
                         + ct_state.dL_x(layout_.template previous<Direction::Z>(idx)));
            auto dE = 0.5
                      * (ct_state.dR_x(idx)
                         + ct_state.dR_x(layout_.template previous<Direction::Z>(idx)));
            auto dB = 0.5
                      * (ct_state.dL_z(idx)
                         + ct_state.dL_z(layout_.template previous<Direction::X>(idx)));
            auto dT = 0.5
                      * (ct_state.dR_z(idx)
                         + ct_state.dR_z(layout_.template previous<Direction::X>(idx)));

            auto [vxW, vxE] = Rec::template reconstruct<Direction::X>(
                ct_state.vt_z(Component::X), idx);
            auto [vzB, vzT] = Rec::template reconstruct<Direction::Z>(
                ct_state.vt_x(Component::Z), idx);
            auto [BzW, BzE]
                = Rec::template reconstruct<Direction::X>(B(Component::Z), idx);
            auto [BxB, BxT]
                = Rec::template reconstruct<Direction::Z>(B(Component::X), idx);

            Ey(idx) = (aW * vxW * BzW + aE * vxE * BzE) - (aB * vzB * BxB + aT * vzT * BxT)
                      - (dE * BzE - dW * BzW) + (dT * BxT - dB * BxB);

            if constexpr (Hall && !IdealOnly)
            {
                auto [jxW, jxE] = Rec::template reconstruct<Direction::X>(
                    ct_state.jt_z(Component::X), idx);
                auto [jzB, jzT] = Rec::template reconstruct<Direction::Z>(
                    ct_state.jt_x(Component::Z), idx);
                auto [rhoW, rhoE]
                    = Rec::template reconstruct<Direction::X>(ct_state.rhot_z, idx);
                auto [rhoB, rhoT]
                    = Rec::template reconstruct<Direction::Z>(ct_state.rhot_x, idx);
                Ey(idx) += -(aW * jxW * BzW / rhoW + aE * jxE * BzE / rhoE)
                           + (aB * jzB * BxB / rhoB + aT * jzT * BxT / rhoT);
            }
        }
    }

    template<typename Rec = Reconstruction_t, bool IdealOnly = false>
    void EzEq_(auto& ct_state, auto& Ez, auto const& B, MeshIndex<dimension> idx) const
    {
        if constexpr (dimension == 1)
        {
            auto [ByL, ByR]
                = Rec::template reconstruct<Direction::X>(B(Component::Y), idx);

            auto FL = ByL * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Y)(idx);
            auto FR = ByR * ct_state.vt_x(Component::X)(idx)
                      - B(Component::X)(idx) * ct_state.vt_x(Component::Y)(idx);

            Ez(idx) = -(ct_state.aL_x(idx) * FL + ct_state.aR_x(idx) * FR)
                      + (ct_state.dR_x(idx) * ByR - ct_state.dL_x(idx) * ByL);

            if constexpr (Hall && !IdealOnly)
            {
                auto invRho  = 1.0 / ct_state.rhot_x(idx);
                auto JxB_z_L = ct_state.jt_x(Component::X)(idx) * ByL
                               - ct_state.jt_x(Component::Y)(idx) * B(Component::X)(idx);
                auto JxB_z_R = ct_state.jt_x(Component::X)(idx) * ByR
                               - ct_state.jt_x(Component::Y)(idx) * B(Component::X)(idx);

                auto HallL = -JxB_z_L * invRho;
                auto HallR = -JxB_z_R * invRho;

                auto F_By_x = ct_state.aL_x(idx) * HallL + ct_state.aR_x(idx) * HallR;

                Ez(idx) += -F_By_x;
            }
        }
        else if constexpr (dimension >= 2)
        {
            auto aW = 0.5
                      * (ct_state.aL_x(idx)
                         + ct_state.aL_x(layout_.template previous<Direction::Y>(idx)));
            auto aE = 0.5
                      * (ct_state.aR_x(idx)
                         + ct_state.aR_x(layout_.template previous<Direction::Y>(idx)));
            auto aS = 0.5
                      * (ct_state.aL_y(idx)
                         + ct_state.aL_y(layout_.template previous<Direction::X>(idx)));
            auto aN = 0.5
                      * (ct_state.aR_y(idx)
                         + ct_state.aR_y(layout_.template previous<Direction::X>(idx)));
            auto dW = 0.5
                      * (ct_state.dL_x(idx)
                         + ct_state.dL_x(layout_.template previous<Direction::Y>(idx)));
            auto dE = 0.5
                      * (ct_state.dR_x(idx)
                         + ct_state.dR_x(layout_.template previous<Direction::Y>(idx)));
            auto dS = 0.5
                      * (ct_state.dL_y(idx)
                         + ct_state.dL_y(layout_.template previous<Direction::X>(idx)));
            auto dN = 0.5
                      * (ct_state.dR_y(idx)
                         + ct_state.dR_y(layout_.template previous<Direction::X>(idx)));

            auto [vyS, vyN] = Rec::template reconstruct<Direction::Y>(
                ct_state.vt_x(Component::Y), idx);
            auto [vxW, vxE] = Rec::template reconstruct<Direction::X>(
                ct_state.vt_y(Component::X), idx);

            auto [BxS, BxN]
                = Rec::template reconstruct<Direction::Y>(B(Component::X), idx);
            auto [ByW, ByE]
                = Rec::template reconstruct<Direction::X>(B(Component::Y), idx);

            Ez(idx) = -(aW * vxW * ByW + aE * vxE * ByE) + (aS * vyS * BxS + aN * vyN * BxN)
                      + (dE * ByE - dW * ByW) - (dN * BxN - dS * BxS);

            if constexpr (Hall && !IdealOnly)
            {
                auto [jyS, jyN] = Rec::template reconstruct<Direction::Y>(
                    ct_state.jt_x(Component::Y), idx);
                auto [jxW, jxE] = Rec::template reconstruct<Direction::X>(
                    ct_state.jt_y(Component::X), idx);

                auto [rhoS, rhoN]
                    = Rec::template reconstruct<Direction::Y>(ct_state.rhot_x, idx);
                auto [rhoW, rhoE]
                    = Rec::template reconstruct<Direction::X>(ct_state.rhot_y, idx);

                Ez(idx) += (aW * jxW * ByW / rhoW + aE * jxE * ByE / rhoE)
                           - (aS * jyS * BxS / rhoS + aN * jyN * BxN / rhoN);
            }
        }
    }

    template<typename Field>
    void resistive_contribution_(Field& E, Field const& J, MeshIndex<Field::dimension> index) const
    {
        E(index) += eta * J(index);
    }

    template<auto component, typename Field, typename VecField>
    void hyperresistive_contribution_(Field& E, Field const& J, VecField const& B, Field const& rho,
                                      MeshIndex<Field::dimension> index) const
    {
        if (hyper_mode == HyperMode::constant)
            return constant_hyperresistive_<component>(E, J, index);
        else if (hyper_mode == HyperMode::spatial)
            return spatial_hyperresistive_<component>(E, J, B, rho, index);
        else
            throw std::runtime_error("Error - Ohm - unknown hyper_mode");
    }

    template<auto component, typename Field>
    void constant_hyperresistive_(Field& E, Field const& J, MeshIndex<Field::dimension> index) const
    {
        E(index) -= nu * layout_.laplacian(J, index);
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

        auto computeHR = [&]<auto BxProj, auto ByProj, auto BzProj, auto rhoProj>() {
            auto const BxOnE = GridLayout::template project<BxProj>(B(Component::X), index);
            auto const ByOnE = GridLayout::template project<ByProj>(B(Component::Y), index);
            auto const BzOnE = GridLayout::template project<BzProj>(B(Component::Z), index);
            auto const nOnE  = GridLayout::template project<rhoProj>(rho, index);
            auto b           = std::sqrt(BxOnE * BxOnE + ByOnE * ByOnE + BzOnE * BzOnE);
            E(index)
                -= nu * layout_.laplacian(J, index) * minMeshSize * minMeshSize * (b / nOnE + 1);
        };

        if constexpr (component == Component::X)
        {
            return computeHR
                .template operator()<GridLayout::BxToEx, GridLayout::ByToEx, GridLayout::BzToEx,
                                     GridLayout::cellCenterToEdgeX>();
        }
        if constexpr (component == Component::Y)
        {
            return computeHR
                .template operator()<GridLayout::BxToEy, GridLayout::ByToEy, GridLayout::BzToEy,
                                     GridLayout::cellCenterToEdgeY>();
        }
        if constexpr (component == Component::Z)
        {
            return computeHR
                .template operator()<GridLayout::BxToEz, GridLayout::ByToEz, GridLayout::BzToEz,
                                     GridLayout::cellCenterToEdgeZ>();
        }
    }


    GridLayout layout_;
};
} // namespace PHARE::core

#endif
