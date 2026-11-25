#ifndef CORE_NUMERICS_RIEMANN_SOLVERS_RUSANOV_HPP
#define CORE_NUMERICS_RIEMANN_SOLVERS_RUSANOV_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"

namespace PHARE::core
{
// gridlayout not used right now as we dont need the inversemeshdir, torm maybe
// do we want to keep supporting whistler waves in the fan ? maybe even then we want a better
// mechanism for the 1/dx
class Rusanov
{
public:
    Rusanov(double const gamma)
        : gamma_{gamma}
    {
    }

    template<auto direction>
    auto solve(auto& uL, auto& uR, auto const& fL, auto const& fR)
    {
        auto const hydro_speed = rusanov_speeds_<direction>(uL, uR);

        uL.to_conservative(gamma_);
        uR.to_conservative(gamma_);

        auto const [Frho, FrhoVx, FrhoVy, FrhoVz, FBx, FBy, FBz, FEtot]
            = rusanov_(uL.as_tuple(), uR.as_tuple(), fL.as_tuple(), fR.as_tuple(), hydro_speed);

        return PerIndex{Frho, {FrhoVx, FrhoVy, FrhoVz}, {FBx, FBy, FBz}, FEtot};
    }

    // the normal component is actually needed for the 1D riemann solver for E in 1D and 2D. We
    // could have an if constexpr on the dimension there, and have an extra template on the
    // direction.
    static auto vector_riemann_averaging(auto const& L, auto const& R)
    {
        return PerIndexVector<double>{0.5 * (L.x + R.x), 0.5 * (L.y + R.y), 0.5 * (L.z + R.z)};
    }


    std::array<double, 4> uct_coefs;
    PerIndexVector<double> vt{std::nan(""), std::nan(""), std::nan("")};

private:
    double const gamma_;

    template<auto direction>
    auto rusanov_speeds_(auto const& uL, auto const& uR)
    {
        auto const BdotBL = uL.B.x * uL.B.x + uL.B.y * uL.B.y + uL.B.z * uL.B.z;
        auto const BdotBR = uR.B.x * uR.B.x + uR.B.y * uR.B.y + uR.B.z * uR.B.z;

        auto compute_speeds = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR,
                                  auto VcompL, auto VcompR, auto BcompL, auto BcompR) {
            auto cfastL = compute_fast_magnetosonic_(gamma_, uL.rho, BcompL, BdotBL, uL.P);
            auto cfastR = compute_fast_magnetosonic_(gamma_, uR.rho, BcompR, BdotBR, uR.P);
            auto S      = std::max(std::abs(VcompL) + cfastL, std::abs(VcompR) + cfastR);

            uct_coefs_(uL, uR, S);

            return S;
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

    auto rusanov_(auto const& uL, auto const& uR, auto const& fL, auto const& fR,
                  auto const S) const
    {
        auto constexpr N_elements = std::tuple_size_v<std::decay_t<decltype(uL)>>;

        return for_N<N_elements, for_N_R_mode::make_tuple>([&](auto i) {
            return (std::get<i>(fL) + std::get<i>(fR)) * 0.5
                   - S * (std::get<i>(uR) - std::get<i>(uL)) * 0.5;
        });
    }

    void uct_coefs_(auto const& uL, auto const& uR, auto const S)
    {
        uct_coefs[0] = 0.5;
        uct_coefs[1] = 0.5;
        uct_coefs[2] = 0.5 * S;
        uct_coefs[3] = 0.5 * S;
        vt           = vector_riemann_averaging(uL.V, uR.V);
    }
};
} // namespace PHARE::core

#endif
