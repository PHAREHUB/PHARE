#ifndef CORE_NUMERICS_RIEMANN_SOLVERS_HLL_HPP
#define CORE_NUMERICS_RIEMANN_SOLVERS_HLL_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"

namespace PHARE::core
{
template<typename GridLayout, bool Hall>
class HLL
{
public:
    HLL(GridLayout const& layout, double const gamma)
        : layout_{layout}
        , gamma_{gamma}
    {
    }

    template<auto direction>
    auto solve(auto& uL, auto& uR, auto const& fL, auto const& fR) const
    {
        auto const [speedsL, speedsR] = hll_speeds_<direction>(uL, uR);

        uL.to_conservative(gamma_);
        uR.to_conservative(gamma_);

        auto const [hydro_speedL, mag_speedL] = speedsL;
        auto const [hydro_speedR, mag_speedR] = speedsR;

        if constexpr (Hall)
        {
            auto split = [](auto const& a) {
                auto hydro = std::make_tuple(a.rho, a.rhoV().x, a.rhoV().y, a.rhoV().z);
                auto mag   = std::make_tuple(a.B.x, a.B.y, a.B.z, a.Etot());
                return std::make_pair(hydro, mag);
            };

            auto [uLhydro, uLmag] = split(uL);
            auto [uRhydro, uRmag] = split(uR);

            auto const [fLhydro, fLmag] = split(fL);
            auto const [fRhydro, fRmag] = split(fR);

            auto [Frho, FrhoVx, FrhoVy, FrhoVz]
                = hll_(uLhydro, uRhydro, fLhydro, fRhydro, hydro_speedL, hydro_speedR);
            auto [FBx, FBy, FBz, FEtot] = hll_(uLmag, uRmag, fLmag, fRmag, mag_speedL, mag_speedR);

            return PerIndex{Frho, {FrhoVx, FrhoVy, FrhoVz}, {FBx, FBy, FBz}, FEtot};
        }
        else
        {
            auto const [Frho, FrhoVx, FrhoVy, FrhoVz, FBx, FBy, FBz, FEtot]
                = hll_(uL.as_tuple(), uR.as_tuple(), fL.as_tuple(), fR.as_tuple(), hydro_speedL,
                       hydro_speedR);

            return PerIndex{Frho, {FrhoVx, FrhoVy, FrhoVz}, {FBx, FBy, FBz}, FEtot};
        }
    }

private:
    GridLayout layout_;
    double const gamma_;

    template<auto direction>
    auto hll_speeds_(auto const& uL, auto const& uR) const
    {
        auto const BdotBL = uL.B.x * uL.B.x + uL.B.y * uL.B.y + uL.B.z * uL.B.z;
        auto const BdotBR = uR.B.x * uR.B.x + uR.B.y * uR.B.y + uR.B.z * uR.B.z;

        auto compute_speeds = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR,
                                  auto VcompL, auto VcompR, auto BcompL, auto BcompR) {
            auto cfastL = compute_fast_magnetosonic_(gamma_, uL.rho, BcompL, BdotBL, uL.P);
            auto cfastR = compute_fast_magnetosonic_(gamma_, uR.rho, BcompR, BdotBR, uR.P);
            auto SL     = std::min(VcompL - cfastL, VcompR - cfastR);
            auto SR     = std::max(VcompL + cfastL, VcompR + cfastR);
            auto SLb    = SL;
            auto SRb    = SR;

            if constexpr (Hall)
            {
                auto cwL = compute_whistler_(layout_.inverseMeshSize(direction), uL.rho, BdotBL);
                auto cwR = compute_whistler_(layout_.inverseMeshSize(direction), uR.rho, BdotBR);
                SLb      = std::min(VcompL - cfastL - cwL, VcompR - cfastR - cwR);
                SRb      = std::max(VcompL + cfastL + cwL, VcompR + cfastR + cwR);
            }

            return std::make_tuple(std::make_tuple(SL, SR), std::make_tuple(SLb, SRb));
        };

        if constexpr (direction == Direction::X)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.x, uR.V.x,
                                  uL.B.x, uR.B.x);
        else if constexpr (direction == Direction::Y)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.y, uR.V.y,
                                  uL.B.y, uR.B.y);
        else if constexpr (direction == Direction::Z)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.z, uR.V.z,
                                  uL.B.z, uR.B.z);
    }

    auto hll_(auto const& uL, auto const& uR, auto const& fL, auto const& fR, auto const& SL,
              auto const& SR) const
    {
        auto constexpr N_elements = std::tuple_size_v<std::decay_t<decltype(uL)>>;

        auto hll = [&](auto const ul, auto const ur, auto const fl, auto const fr) {
            if (SL > 0.0)
                return fl;
            else if (SR < 0.0)
                return fr;
            else
                return (SR * fl - SL * fr + SL * SR * (ur - ul)) / (SR - SL);
        };

        return for_N<N_elements, for_N_R_mode::make_tuple>([&](auto i) {
            return hll(std::get<i>(uL), std::get<i>(uR), std::get<i>(fL), std::get<i>(fR));
        });
    }
};
} // namespace PHARE::core

#endif
