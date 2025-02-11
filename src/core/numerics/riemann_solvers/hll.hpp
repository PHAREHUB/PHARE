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
    PerIndex solve(PerIndex& uL, PerIndex& uR, PerIndex const& fL, PerIndex const& fR) const
    {
        auto const [sl, sr] = hll_speeds_<direction>(uL, uR);

        uL.toConservative(gamma_);
        uR.toConservative(gamma_);

        if constexpr (Hall)
        {
            auto const [SL, SLb] = sl;
            auto const [SR, SRb] = sr;

            auto ul  = std::make_tuple(uL.rho, uL.Vx, uL.Vy, uL.Vz);
            auto ulb = std::make_tuple(uL.Bx, uL.By, uL.Bz, uL.P);
            auto ur  = std::make_tuple(uR.rho, uR.Vx, uR.Vy, uR.Vz);
            auto urb = std::make_tuple(uR.Bx, uR.By, uR.Bz, uR.P);

            auto fl  = std::make_tuple(fL.rho, fL.Vx, fL.Vy, fL.Vz);
            auto flb = std::make_tuple(fL.Bx, fL.By, fL.Bz, fL.P);
            auto fr  = std::make_tuple(fR.rho, fR.Vx, fR.Vy, fR.Vz);
            auto frb = std::make_tuple(fR.Bx, fR.By, fR.Bz, fR.P);

            auto [Frho, FrhoVx, FrhoVy, FrhoVz] = hll_(ul, ur, fl, fr, SL, SR);
            auto [FBx, FBy, FBz, FEtot]         = hll_(ulb, urb, flb, frb, SLb, SRb);

            return PerIndex{Frho, FrhoVx, FrhoVy, FrhoVz, FBx, FBy, FBz, FEtot};
        }
        else
        {
            auto const [Frho, FrhoVx, FrhoVy, FrhoVz, FBx, FBy, FBz, FEtot]
                = hll_(uL(), uR(), fL(), fR(), std::get<0>(sl), std::get<0>(sr));

            return PerIndex{Frho, FrhoVx, FrhoVy, FrhoVz, FBx, FBy, FBz, FEtot};
        }
    }

private:
    GridLayout layout_;
    double const gamma_;

    template<auto direction>
    auto hll_speeds_(auto const& uL, auto const& uR) const
    {
        auto const BdotBL = uL.Bx * uL.Bx + uL.By * uL.By + uL.Bz * uL.Bz;
        auto const BdotBR = uR.Bx * uR.Bx + uR.By * uR.By + uR.Bz * uR.Bz;

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
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.Vx, uR.Vx, uL.Bx,
                                  uR.Bx);
        else if constexpr (direction == Direction::Y)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.Vy, uR.Vy, uL.By,
                                  uR.By);
        else if constexpr (direction == Direction::Z)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.Vz, uR.Vz, uL.Bz,
                                  uR.Bz);
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
