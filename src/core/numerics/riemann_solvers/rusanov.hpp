#ifndef CORE_NUMERICS_RIEMANN_SOLVERS_RUSANOV_HPP
#define CORE_NUMERICS_RIEMANN_SOLVERS_RUSANOV_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"

namespace PHARE::core
{
// gridlayout not used right now as we dont need the inversemeshdir, torm maybe
// do we want to keep supporting whistler waves in the fan ? maybe even then we want a better
// mechanism for the 1/dx
template<bool Hall>
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
            = rusanov_(uL.as_reduced_conservative_tuple(), uR.as_reduced_conservative_tuple(),
                       fL.as_tuple(), fR.as_tuple(), hydro_speed);

        return PerIndex{Frho, {FrhoVx, FrhoVy, FrhoVz}, {FBx, FBy, FBz}, FEtot};
    }

    template<auto direction>
    auto solve(auto& uL, auto& uR, auto const& fL, auto const& fR, auto const& jL, auto const& jR)
    {
        auto const speeds = rusanov_speeds_<direction>(uL, uR, jL, jR);

        uL.to_conservative(gamma_);
        uR.to_conservative(gamma_);

        auto const [hydro_speed, mag_speed] = speeds;

        auto split_state = [](auto const& a) {
            auto const reduced = a.as_reduced_conservative_tuple();
            auto hydro = std::make_tuple(std::get<0>(reduced), std::get<1>(reduced),
                                         std::get<2>(reduced), std::get<3>(reduced));
            auto mag   = std::make_tuple(std::get<4>(reduced), std::get<5>(reduced),
                                         std::get<6>(reduced), std::get<7>(reduced));
            return std::make_pair(hydro, mag);
        };
        auto split_flux = [](auto const& a) {
            auto const raw = a.as_tuple();
            auto hydro = std::make_tuple(std::get<0>(raw), std::get<1>(raw), std::get<2>(raw),
                                         std::get<3>(raw));
            auto mag
                = std::make_tuple(std::get<4>(raw), std::get<5>(raw), std::get<6>(raw), std::get<7>(raw));
            return std::make_pair(hydro, mag);
        };

        auto [uLhydro, uLmag] = split_state(uL);
        auto [uRhydro, uRmag] = split_state(uR);

        auto const [fLhydro, fLmag] = split_flux(fL);
        auto const [fRhydro, fRmag] = split_flux(fR);

        auto [Frho, FrhoVx, FrhoVy, FrhoVz]
            = rusanov_(uLhydro, uRhydro, fLhydro, fRhydro, hydro_speed);
        auto [FBx, FBy, FBz, FEtot] = rusanov_(uLmag, uRmag, fLmag, fRmag, mag_speed);

        return PerIndex{Frho, {FrhoVx, FrhoVy, FrhoVz}, {FBx, FBy, FBz}, FEtot};
    }

    static auto riemann_averaging(auto const& L, auto const& R) { return 0.5 * (L + R); }

    // the normal component is actually needed for the 1D riemann solver for E in 1D and 2D. We
    // could have an if constexpr on the dimension there, and have an extra template on the
    // direction.
    static auto vector_riemann_averaging(auto const& L, auto const& R)
    {
        return PerIndexVector<double>{0.5 * (L.x + R.x), 0.5 * (L.y + R.y), 0.5 * (L.z + R.z)};
    }


    std::array<double, 4> uct_coefs;
    PerIndexVector<double> vt{std::nan(""), std::nan(""), std::nan("")};

    PerIndexVector<double> jt{std::nan(""), std::nan(""), std::nan("")};
    double rhot{std::nan("")};

private:
    double const gamma_;

    template<auto direction>
    auto rusanov_speeds_(auto const& uL, auto const& uR)
    {
        auto const BtL    = uL.totalB();
        auto const BtR    = uR.totalB();
        auto const BdotBL = BtL.x * BtL.x + BtL.y * BtL.y + BtL.z * BtL.z;
        auto const BdotBR = BtR.x * BtR.x + BtR.y * BtR.y + BtR.z * BtR.z;

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
                                  BtL.x, BtR.x);
        else if constexpr (direction == Direction::Y)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.y, uR.V.y,
                                  BtL.y, BtR.y);
        else if constexpr (direction == Direction::Z)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.z, uR.V.z,
                                  BtL.z, BtR.z);
    }

    template<auto direction>
    auto rusanov_speeds_(auto const& uL, auto const& uR, auto const& jL, auto const& jR)
    {
        auto const BtL    = uL.totalB();
        auto const BtR    = uR.totalB();
        auto const BdotBL = BtL.x * BtL.x + BtL.y * BtL.y + BtL.z * BtL.z;
        auto const BdotBR = BtR.x * BtR.x + BtR.y * BtR.y + BtR.z * BtR.z;

        auto compute_speeds = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR,
                                  auto VcompL, auto VcompR, auto BcompL, auto BcompR) {
            auto cfastL = compute_fast_magnetosonic_(gamma_, uL.rho, BcompL, BdotBL, uL.P);
            auto cfastR = compute_fast_magnetosonic_(gamma_, uR.rho, BcompR, BdotBR, uR.P);

            auto S = std::max(std::abs(VcompL) + cfastL, std::abs(VcompR) + cfastR);

            // no-whistler test
            // might need to do something for the layout if we want to use this again
            auto cwL = 0.; // compute_whistler_(layout_.inverseMeshSize(direction), uL.rho, BdotBL);
            auto cwR = 0.; // compute_whistler_(layout_.inverseMeshSize(direction), uR.rho, BdotBR);

            auto Sb = std::max(std::abs(VcompL) + cfastL + cwL, std::abs(VcompR) + cfastR + cwR);

            uct_coefs_(uL, uR, jL, jR, Sb);

            return std::make_pair(S, Sb);
        };

        if constexpr (direction == Direction::X)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.x, uR.V.x,
                                  BtL.x, BtR.x);
        else if constexpr (direction == Direction::Y)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.y, uR.V.y,
                                  BtL.y, BtR.y);
        else if constexpr (direction == Direction::Z)
            return compute_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.z, uR.V.z,
                                  BtL.z, BtR.z);
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

    void uct_coefs_(auto const& uL, auto const& uR, auto const& jL, auto const& jR, auto const S)
    {
        uct_coefs_(uL, uR, S);

        jt = vector_riemann_averaging(jL, jR);

        rhot = riemann_averaging(uL.rho, uR.rho);
    }
};
} // namespace PHARE::core

#endif
