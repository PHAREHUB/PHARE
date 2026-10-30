#ifndef CORE_NUMERICS_RIEMANN_SOLVERS_HLLD_HPP
#define CORE_NUMERICS_RIEMANN_SOLVERS_HLLD_HPP

#include "core/def.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/riemann_solvers/mhd_speeds.hpp"

#include <cstdlib>
#include <type_traits>

namespace PHARE::core
{


template<bool Hall>
class HLLD
{
public:
    HLLD(double const gamma)
        : gamma_{gamma}
    {
    }

    template<auto direction>
    PHARE_ALWAYS_INLINE auto solve(auto& uL, auto& uR, auto const& fL, auto const& fR)
    {
        auto hlld_speeds = hlld_speeds_<direction>(uL, uR);

        return fluxes_<direction>(uL, uR, fL, fR, hlld_speeds);
    }

    template<auto direction>
    PHARE_ALWAYS_INLINE auto solve(auto& uL, auto& uR, auto const& fL, auto const& fR,
                                   auto const& jL, auto const& jR)
    {
        auto hlld_speeds = hlld_speeds_<direction>(uL, uR, jL, jR);

        return fluxes_<direction>(uL, uR, fL, fR, hlld_speeds);
    }

    // using HLL averages here following PLUTO and Idefix implementation
    auto riemann_averaging(auto const& L, auto const& R) const
    {
        auto const inv = 1.0 / (SR_ - SL_);
        return (SR_ * L - SL_ * R) * inv;
    }

    auto vector_riemann_averaging(auto const& L, auto const& R) const
    {
        auto const inv = 1.0 / (SR_ - SL_);
        return PerIndexVector<double>{(SR_ * L.x - SL_ * R.x) * inv, (SR_ * L.y - SL_ * R.y) * inv,
                                      (SR_ * L.z - SL_ * R.z) * inv};
    }

    std::array<double, 4> uct_coefs;
    PerIndexVector<double> vt{std::nan(""), std::nan(""), std::nan("")};

    PerIndexVector<double> jt{std::nan(""), std::nan(""), std::nan("")};
    double rhot{std::nan("")};

private:
    double const gamma_;

    // the star states reach the flux through a continuation rather than a returned tuple, which
    // GCC keeps in memory and which then blocks the face loop's vectorisation
    template<auto direction>
    PHARE_ALWAYS_INLINE auto fluxes_(auto& uL, auto& uR, auto const& fL, auto const& fR,
                                     auto& hlld_speeds)
    {
        return hlld_intermediate_states_<direction>(
            uL, uR, fL, fR, hlld_speeds,
            [&](auto const& uL_s, auto const& uL_ss, auto const& uR_ss, auto const& uR_s)
                PHARE_LAMBDA_ALWAYS_INLINE {
                    uL.to_conservative(gamma_);
                    uR.to_conservative(gamma_);

                    return hlld_(uL, uL_s, uL_ss, uR_ss, uR_s, uR, fL, fR, hlld_speeds);
                });
    }

    // these are used for the riemann averagings that are always done after speed computation per
    // index. This is needed for the save transverse magnetic field in godunov fluxes for any
    // riemann solver, but maybe not the best interface.
    double SL_;
    double SR_;

    template<auto direction>
    PHARE_ALWAYS_INLINE auto hlld_speeds_(auto const& uL, auto const& uR)
    {
        auto const BdotBL = uL.B.x * uL.B.x + uL.B.y * uL.B.y + uL.B.z * uL.B.z;
        auto const BdotBR = uR.B.x * uR.B.x + uR.B.y * uR.B.y + uR.B.z * uR.B.z;

        auto compute_hll_speeds
            = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR, auto VcompL,
                  auto VcompR, auto BcompL, auto BcompR) PHARE_LAMBDA_ALWAYS_INLINE {
                  auto cfastL = compute_fast_magnetosonic_(gamma_, uL.rho, BcompL, BdotBL, uL.P);
                  auto cfastR = compute_fast_magnetosonic_(gamma_, uR.rho, BcompR, BdotBR, uR.P);
                  // two-argument min/max: the initializer_list overloads keep their operands in
                  // memory, which stops clang vectorising the face loop
                  auto SL = std::min(VcompL - cfastL, VcompR - cfastR);
                  auto SR = std::max(VcompL + cfastL, VcompR + cfastR);

                  return std::make_tuple(SL, SR);
              };

        auto compute_hlld_speeds
            = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR, auto VcompL,
                  auto VcompR, auto BcompL, auto BcompR) PHARE_LAMBDA_ALWAYS_INLINE {
                  auto [SL, SR] = compute_hll_speeds(rhoL, rhoR, PL, PR, BdotBL, BdotBR, VcompL,
                                                     VcompR, BcompL, BcompR);

                  auto PtL = PL + 0.5 * BdotBL;
                  auto PtR = PR + 0.5 * BdotBR;

                  auto const Bn = BcompL; // is the same on both sides
                  // auto const Bn = (SR * BcompR - SL * BcompL) / (SR - SL);

                  auto SM_numerator
                      = rhoR * VcompR * (SR - VcompR) - rhoL * VcompL * (SL - VcompL) - PtR + PtL;
                  auto SM_denominator = rhoR * (SR - VcompR) - rhoL * (SL - VcompL);
                  auto SM             = SM_numerator / SM_denominator;

                  auto const rhoL_s = rhoL * (SL - VcompL) / (SL - SM);
                  auto const rhoR_s = rhoR * (SR - VcompR) / (SR - SM);

                  auto const SL_s = SM - std::abs(Bn) / std::sqrt(rhoL_s);
                  auto const SR_s = SM + std::abs(Bn) / std::sqrt(rhoR_s);

                  auto hlld_speeds = std::make_tuple(SL, SL_s, SM, SR_s, SR);

                  uct_coefs_<direction>(uL, uR, hlld_speeds);


                  return hlld_speeds;
              };


        if constexpr (direction == Direction::X)
            return compute_hlld_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.x, uR.V.x,
                                       uL.B.x, uR.B.x);
        else if constexpr (direction == Direction::Y)
            return compute_hlld_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.y, uR.V.y,
                                       uL.B.y, uR.B.y);
        else if constexpr (direction == Direction::Z)
            return compute_hlld_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.z, uR.V.z,
                                       uL.B.z, uR.B.z);
    }

    // need some optimization, the star states should only be computed if we are in a star region
    template<auto direction>
    PHARE_ALWAYS_INLINE auto hlld_intermediate_states_(auto const& uL, auto const& uR,
                                                       auto const& fL, auto const& fR,
                                                       auto& hlld_speeds, auto&& then) const
    {
        auto& [SL, SL_s, SM, SR_s, SR] = hlld_speeds;

        auto const etotL
            = eosPToEtot(gamma_, uL.rho, uL.V.x, uL.V.y, uL.V.z, uL.B.x, uL.B.y, uL.B.z, uL.P);

        auto const etotR
            = eosPToEtot(gamma_, uR.rho, uR.V.x, uR.V.y, uR.V.z, uR.B.x, uR.B.y, uR.B.z, uR.P);

        auto const p_tot_L = uL.P + 0.5 * (uL.B.x * uL.B.x + uL.B.y * uL.B.y + uL.B.z * uL.B.z);
        auto const p_tot_R = uR.P + 0.5 * (uR.B.x * uR.B.x + uR.B.y * uR.B.y + uR.B.z * uR.B.z);


        auto constexpr transverse = [&]() PHARE_LAMBDA_ALWAYS_INLINE {
            if constexpr (direction == Direction::X)
                return std::array{Direction::Y, Direction::Z};
            else if constexpr (direction == Direction::Y)
                return std::array{Direction::X, Direction::Z};
            else if constexpr (direction == Direction::Z)
                return std::array{Direction::X, Direction::Y};
        }();

        auto sgn = [](auto const x) PHARE_LAMBDA_ALWAYS_INLINE { return (x > 0.0) - (x < 0.0); };

        auto const inv_SRmSL = 1.0 / (SR - SL);

        auto compute_tranverse_magnetic_s_hllc = [&](auto const tdir) PHARE_LAMBDA_ALWAYS_INLINE {
            auto const Bhll
                = (SR * uR.B(tdir) - SL * uL.B(tdir) + fL.B(tdir) - fR.B(tdir)) * inv_SRmSL;

            return Bhll;
        };

        auto hlld_intermediate_states = [&](auto const fl,
                                            auto const fr) PHARE_LAMBDA_ALWAYS_INLINE {
            auto const rhoL_s = uL.rho * (SL - uL.V(direction)) / (SL - SM);
            auto const rhoR_s = uR.rho * (SR - uR.V(direction)) / (SR - SM);

            auto const vn = SM;
            // this should probably not be reconstructed in the normal direction as B is already
            // face centered there
            auto const Bn = uL.B(direction); // should be the same on both sides
            // auto const Bn = (SR * uR.B(direction) - SL * uL.B(direction)) / (SR - SL);


            // the factor scaling both transverse components of B on one side
            auto transverse_magnetic_factor
                = [&](auto const& u, auto const S) PHARE_LAMBDA_ALWAYS_INLINE {
                      return (u.rho * (S - u.V(direction)) * (S - u.V(direction)) - Bn * Bn)
                             / (u.rho * (S - u.V(direction)) * (S - SM) - Bn * Bn);
                  };
            auto const btL = transverse_magnetic_factor(uL, SL);
            auto const btR = transverse_magnetic_factor(uR, SR);

            auto Bt0L_s = uL.B(transverse[0]) * btL;
            auto Bt1L_s = uL.B(transverse[1]) * btL;
            auto Bt0R_s = uR.B(transverse[0]) * btR;
            auto Bt1R_s = uR.B(transverse[1]) * btR;

            // | rather than ||: both comparisons always run, so neither is a branch
            bool const hllc_fallback
                = ((SL_s - SL) < 1.0e-4 * (SM - SL)) | ((SR_s - SR) > -1.0e-4 * (SR - SM));

            // Fallback to HLLC as in idefix/PLUTO (Maybe Bn should also use HLL averaging?).
            // Both branches are computed and selected, so the face loop has no control flow and
            // vectorises.
            auto const Bn_hllc  = compute_tranverse_magnetic_s_hllc(direction);
            auto const Bt0_hllc = compute_tranverse_magnetic_s_hllc(transverse[0]);
            auto const Bt1_hllc = compute_tranverse_magnetic_s_hllc(transverse[1]);

            auto const Bn_s = hllc_fallback ? Bn_hllc : Bn;
            Bt0L_s          = hllc_fallback ? Bt0_hllc : Bt0L_s;
            Bt0R_s          = hllc_fallback ? Bt0_hllc : Bt0R_s;
            Bt1L_s          = hllc_fallback ? Bt1_hllc : Bt1L_s;
            Bt1R_s          = hllc_fallback ? Bt1_hllc : Bt1R_s;
            SL_s            = hllc_fallback ? SM : SL_s;
            SR_s            = hllc_fallback ? SM : SR_s;

            // auto compute_tranverse_velocity_s = [&](auto const u, auto const S, auto const tdir)
            // {
            //     auto const vt_s = u.V(tdir)
            //                       - u.B(tdir) * Bn * (SM - u.V(direction))
            //                             / (u.rho * (S - u.V(direction)) * (S - SM) - Bn * Bn);
            //
            //     return vt_s;
            // };
            //
            // auto const vt0L_s = compute_tranverse_velocity_s(uL, SL, transverse[0]);
            // auto const vt1L_s = compute_tranverse_velocity_s(uL, SL, transverse[1]);
            // auto const vt0R_s = compute_tranverse_velocity_s(uR, SR, transverse[0]);
            // auto const vt1R_s = compute_tranverse_velocity_s(uR, SR, transverse[1]);

            auto const BnL = Bn / (uL.rho * (SL - uL.V(direction)));
            auto const BnR = Bn / (uR.rho * (SR - uR.V(direction)));

            auto const vt0L_s = uL.V(transverse[0]) - (Bt0L_s - uL.B(transverse[0])) * BnL;
            auto const vt1L_s = uL.V(transverse[1]) - (Bt1L_s - uL.B(transverse[1])) * BnL;
            auto const vt0R_s = uR.V(transverse[0]) - (Bt0R_s - uR.B(transverse[0])) * BnR;
            auto const vt1R_s = uR.V(transverse[1]) - (Bt1R_s - uR.B(transverse[1])) * BnR;


            auto const p_tot_s
                = ((SR - uR.V(direction)) * uR.rho * p_tot_L
                   - (SL - uL.V(direction)) * uL.rho * p_tot_R
                   + uL.rho * uR.rho * (SL - uL.V(direction)) * (SR - uR.V(direction))
                         * (uR.V(direction) - uL.V(direction)))
                  / (uR.rho * (SR - uR.V(direction)) - uL.rho * (SL - uL.V(direction)));

            auto const EtotL_s
                = ((SL - uL.V(direction)) * etotL - p_tot_L * uL.V(direction) + p_tot_s * SM
                   + Bn
                         * (Bn * uL.V(direction) + uL.B(transverse[0]) * uL.V(transverse[0])
                            + uL.B(transverse[1]) * uL.V(transverse[1])
                            - (Bn * vn + vt0L_s * Bt0L_s + vt1L_s * Bt1L_s)))
                  / (SL - SM);
            auto const EtotR_s
                = ((SR - uR.V(direction)) * etotR - p_tot_R * uR.V(direction) + p_tot_s * SM
                   + Bn
                         * (Bn * uR.V(direction) + uR.B(transverse[0]) * uR.V(transverse[0])
                            + uR.B(transverse[1]) * uR.V(transverse[1])
                            - (Bn * vn + vt0R_s * Bt0R_s + vt1R_s * Bt1R_s)))
                  / (SR - SM);


            auto const sqrtL   = std::sqrt(rhoL_s);
            auto const sqrtR   = std::sqrt(rhoR_s);
            auto const sqrtLR  = sqrtL * sqrtR;
            auto const inv_sum = 1.0 / (sqrtL + sqrtR);

            auto const vt0_ss
                = (sqrtL * vt0L_s + sqrtR * vt0R_s + (Bt0R_s - Bt0L_s) * sgn(Bn)) * inv_sum;
            auto const vt1_ss
                = (sqrtL * vt1L_s + sqrtR * vt1R_s + (Bt1R_s - Bt1L_s) * sgn(Bn)) * inv_sum;

            auto const Bt0_ss
                = (sqrtL * Bt0R_s + sqrtR * Bt0L_s + sqrtLR * (vt0R_s - vt0L_s) * sgn(Bn))
                  * inv_sum;
            auto const Bt1_ss
                = (sqrtL * Bt1R_s + sqrtR * Bt1L_s + sqrtLR * (vt1R_s - vt1L_s) * sgn(Bn))
                  * inv_sum;

            auto const EtotL_ss = EtotL_s
                                  - sqrtL
                                        * (vn * Bn + vt0L_s * Bt0L_s + vt1L_s * Bt1L_s
                                           - (vn * Bn + vt0_ss * Bt0_ss + vt1_ss * Bt1_ss))
                                        * sgn(Bn);
            auto const EtotR_ss = EtotR_s
                                  + sqrtR
                                        * (vn * Bn + vt0R_s * Bt0R_s + vt1R_s * Bt1R_s
                                           - (vn * Bn + vt0_ss * Bt0_ss + vt1_ss * Bt1_ss))
                                        * sgn(Bn);

            PerIndexVector<std::decay_t<decltype(rhoL_s)>> rhoVL_s, rhoVR_s, rhoVL_ss, rhoVR_ss;
            PerIndexVector<std::decay_t<decltype(rhoL_s)>> BL_s, BR_s, BL_ss, BR_ss;

            rhoVL_s(direction)      = rhoL_s * vn;
            rhoVL_s(transverse[0])  = rhoL_s * vt0L_s;
            rhoVL_s(transverse[1])  = rhoL_s * vt1L_s;
            rhoVR_s(direction)      = rhoR_s * vn;
            rhoVR_s(transverse[0])  = rhoR_s * vt0R_s;
            rhoVR_s(transverse[1])  = rhoR_s * vt1R_s;
            rhoVL_ss(direction)     = rhoL_s * vn;
            rhoVL_ss(transverse[0]) = rhoL_s * vt0_ss;
            rhoVL_ss(transverse[1]) = rhoL_s * vt1_ss;
            rhoVR_ss(direction)     = rhoR_s * vn;
            rhoVR_ss(transverse[0]) = rhoR_s * vt0_ss;
            rhoVR_ss(transverse[1]) = rhoR_s * vt1_ss;

            BL_s(direction)      = Bn_s;
            BL_s(transverse[0])  = Bt0L_s;
            BL_s(transverse[1])  = Bt1L_s;
            BR_s(direction)      = Bn_s;
            BR_s(transverse[0])  = Bt0R_s;
            BR_s(transverse[1])  = Bt1R_s;
            BL_ss(direction)     = Bn;
            BL_ss(transverse[0]) = Bt0_ss;
            BL_ss(transverse[1]) = Bt1_ss;
            BR_ss(direction)     = Bn;
            BR_ss(transverse[0]) = Bt0_ss;
            BR_ss(transverse[1]) = Bt1_ss;

            auto const uL_s  = PerIndex{rhoL_s, rhoVL_s, BL_s, EtotL_s};
            auto const uR_s  = PerIndex{rhoR_s, rhoVR_s, BR_s, EtotR_s};
            auto const uL_ss = PerIndex{rhoL_s, rhoVL_ss, BL_ss, EtotL_ss};
            auto const uR_ss = PerIndex{rhoR_s, rhoVR_ss, BR_ss, EtotR_ss};

            return then(uL_s, uL_ss, uR_ss, uR_s);
        };

        return hlld_intermediate_states(fL, fR);
    }

    // template<auto direction>
    // auto hlld_intermediate_states_(auto const& uL, auto const& uR, auto const& fL, auto const&
    // fR,
    //                                auto& hlld_speeds) const
    // {
    //     auto& [SL, SL_s, SM, SR_s, SR] = hlld_speeds;
    //
    //     auto const etotL
    //         = eosPToEtot(gamma_, uL.rho, uL.V.x, uL.V.y, uL.V.z, uL.B.x, uL.B.y, uL.B.z, uL.P);
    //
    //     auto const etotR
    //         = eosPToEtot(gamma_, uR.rho, uR.V.x, uR.V.y, uR.V.z, uR.B.x, uR.B.y, uR.B.z, uR.P);
    //
    //     auto const p_tot_L = uL.P + 0.5 * (uL.B.x * uL.B.x + uL.B.y * uL.B.y + uL.B.z * uL.B.z);
    //     auto const p_tot_R = uR.P + 0.5 * (uR.B.x * uR.B.x + uR.B.y * uR.B.y + uR.B.z * uR.B.z);
    //
    //
    //     auto constexpr transverse = [&]() {
    //         if constexpr (direction == Direction::X)
    //             return std::array{Direction::Y, Direction::Z};
    //         else if constexpr (direction == Direction::Y)
    //             return std::array{Direction::X, Direction::Z};
    //         else if constexpr (direction == Direction::Z)
    //             return std::array{Direction::X, Direction::Y};
    //     }();
    //
    //     auto hlld_intermediate_states = [&](auto const fl, auto const fr) {
    //         auto const sLmV = SL - uL.V(direction);
    //         auto const sRmV = SR - uR.V(direction);
    //
    //         auto const sLmM = SL - SM;
    //         auto const sRmM = SR - SM;
    //
    //         auto const rhoL_s = uL.rho * sLmV / sLmM;
    //         auto const rhoR_s = uR.rho * sRmV / sRmM;
    //
    //         auto const Bn = uL.B(direction); // should be the same on both sides
    //
    //         auto const sgnBn = Bn > 0. ? 1. : -1.;
    //
    //
    //         auto const vn = SM;
    //
    //         auto const p_tot_L_s = p_tot_L + uL.rho * sLmV * (SM - uL.V(direction));
    //         auto const p_tot_R_s = p_tot_R + uR.rho * sRmV * (SM - uR.V(direction));
    //
    //         auto const p_tot_s = 0.5 * (p_tot_L_s + p_tot_R_s);
    //
    //         using Float = std::decay_t<decltype(rhoL_s)>;
    //         PerIndexVector<Float> rhoVL_s, BL_s;
    //
    //         rhoVL_s(direction) = rhoL_s * vn;
    //         BL_s(direction)    = Bn;
    //         if (std::abs(uL.rho * sLmV * sLmM - Bn * Bn) < 1e-4 * p_tot_s)
    //         {
    //             // degenerate case
    //             rhoVL_s(transverse[0]) = rhoL_s * uL.V(transverse[0]);
    //             rhoVL_s(transverse[1]) = rhoL_s * uL.V(transverse[1]);
    //
    //             BL_s(transverse[0]) = uL.B(transverse[0]);
    //             BL_s(transverse[1]) = uL.B(transverse[1]);
    //         }
    //         else
    //         {
    //             // normal case
    //             auto tmp               = Bn * (sLmV - sLmM) / (uL.rho * sLmV * sLmM - Bn * Bn);
    //             rhoVL_s(transverse[0]) = rhoL_s * (uL.V(transverse[0]) - uL.B(transverse[0]) *
    //             tmp); rhoVL_s(transverse[1]) = rhoL_s * (uL.V(transverse[1]) -
    //             uL.B(transverse[1]) * tmp);
    //
    //             tmp = (uL.rho * sLmV * sLmV - Bn * Bn) / (uL.rho * sLmV * sLmM - Bn * Bn);
    //             BL_s(transverse[0]) = tmp * uL.B(transverse[0]);
    //             BL_s(transverse[1]) = tmp * uL.B(transverse[1]);
    //         }
    //
    //         auto const vdBL_s
    //             = (rhoVL_s(direction) * Bn + rhoVL_s(transverse[0]) * BL_s(transverse[0])
    //                + rhoVL_s(transverse[1]) * BL_s(transverse[1]))
    //               / rhoL_s;
    //
    //         auto const EtotL_s
    //             = (sLmV * etotL - p_tot_L * uL.V(direction) + p_tot_s * SM
    //                + Bn
    //                      * (Bn * uL.V(direction) + uL.B(transverse[0]) * uL.V(transverse[0])
    //                         + uL.B(transverse[1]) * uL.V(transverse[1]) - vdBL_s))
    //               / sLmM;
    //
    //
    //         PerIndexVector<Float> rhoVR_s, BR_s;
    //
    //         rhoVR_s(direction) = rhoR_s * vn;
    //         BR_s(direction)    = Bn;
    //         if (std::abs(uR.rho * sRmV * sRmM - Bn * Bn) < 1e-4 * p_tot_s)
    //         {
    //             // degenerate case
    //             rhoVR_s(transverse[0]) = rhoR_s * uR.V(transverse[0]);
    //             rhoVR_s(transverse[1]) = rhoR_s * uR.V(transverse[1]);
    //
    //             BR_s(transverse[0]) = uR.B(transverse[0]);
    //             BR_s(transverse[1]) = uR.B(transverse[1]);
    //         }
    //         else
    //         {
    //             // normal case
    //             auto tmp               = Bn * (sRmV - sRmM) / (uR.rho * sRmV * sRmM - Bn * Bn);
    //             rhoVR_s(transverse[0]) = rhoR_s * (uR.V(transverse[0]) - uR.B(transverse[0]) *
    //             tmp); rhoVR_s(transverse[1]) = rhoR_s * (uR.V(transverse[1]) -
    //             uR.B(transverse[1]) * tmp);
    //
    //             tmp = (uR.rho * sRmV * sRmV - Bn * Bn) / (uR.rho * sRmV * sRmM - Bn * Bn);
    //             BR_s(transverse[0]) = tmp * uR.B(transverse[0]);
    //             BR_s(transverse[1]) = tmp * uR.B(transverse[1]);
    //         }
    //
    //         auto const vdBR_s
    //             = (rhoVR_s(direction) * Bn + rhoVR_s(transverse[0]) * BR_s(transverse[0])
    //                + rhoVR_s(transverse[1]) * BR_s(transverse[1]))
    //               / rhoR_s;
    //
    //         auto const EtotR_s
    //             = (sRmV * etotR - p_tot_R * uR.V(direction) + p_tot_s * SM
    //                + Bn
    //                      * (Bn * uR.V(direction) + uR.B(transverse[0]) * uR.V(transverse[0])
    //                         + uR.B(transverse[1]) * uR.V(transverse[1]) - vdBR_s))
    //               / sRmM;
    //
    //
    //         PerIndexVector<Float> rhoVL_ss, rhoVR_ss, BL_ss, BR_ss;
    //
    //         rhoVL_ss(direction) = rhoL_s * vn;
    //         rhoVR_ss(direction) = rhoR_s * vn;
    //         BL_ss(direction) = BR_ss(direction) = Bn;
    //
    //         auto const sqrt_rhoL_s      = std::sqrt(rhoL_s);
    //         auto const sqrt_rhoR_s      = std::sqrt(rhoR_s);
    //         auto const inv_sqrt_rho_sum = 1.0 / (sqrt_rhoL_s + sqrt_rhoR_s);
    //
    //         auto tmp = inv_sqrt_rho_sum
    //                    * (sqrt_rhoL_s * (rhoVL_s(transverse[0]) / rhoL_s)
    //                       + sqrt_rhoR_s * (rhoVR_s(transverse[0]) / rhoR_s));
    //
    //         rhoVL_ss(transverse[0]) = rhoL_s * tmp;
    //         rhoVR_ss(transverse[0]) = rhoR_s * tmp;
    //
    //         tmp = inv_sqrt_rho_sum
    //               * (sqrt_rhoL_s * (rhoVL_s(transverse[1]) / rhoL_s)
    //                  + sqrt_rhoR_s * (rhoVR_s(transverse[1]) / rhoR_s));
    //
    //         rhoVL_ss(transverse[1]) = rhoL_s * tmp;
    //         rhoVR_ss(transverse[1]) = rhoR_s * tmp;
    //
    //         tmp = inv_sqrt_rho_sum
    //               * (sqrt_rhoL_s * BL_s(transverse[0]) + sqrt_rhoR_s * BR_s(transverse[0])
    //                  + sgnBn * sqrt_rhoL_s * sqrt_rhoR_s
    //                        * (rhoVR_s(transverse[0]) / rhoR_s - rhoVL_s(transverse[0]) /
    //                        rhoL_s));
    //
    //         BL_ss(transverse[0]) = BR_ss(transverse[0]) = tmp;
    //
    //         tmp = inv_sqrt_rho_sum
    //               * (sqrt_rhoL_s * BL_s(transverse[1]) + sqrt_rhoR_s * BR_s(transverse[1])
    //                  + sgnBn * sqrt_rhoL_s * sqrt_rhoR_s
    //                        * (rhoVR_s(transverse[1]) / rhoR_s - rhoVL_s(transverse[1]) /
    //                        rhoL_s));
    //
    //         BL_ss(transverse[1]) = BR_ss(transverse[1]) = tmp;
    //
    //         tmp = SM * Bn + (rhoVL_ss(transverse[0]) * BL_ss(transverse[0])) / rhoL_s
    //               + (rhoVL_ss(transverse[1]) * BL_ss(transverse[1])) / rhoL_s;
    //
    //         auto const EtotL_ss = EtotL_s - sqrt_rhoL_s * sgnBn * (vdBL_s - tmp);
    //         auto const EtotR_ss = EtotR_s + sqrt_rhoR_s * sgnBn * (vdBR_s - tmp);
    //
    //         auto const uL_s  = PerIndex{rhoL_s, rhoVL_s, BL_s, EtotL_s};
    //         auto const uR_s  = PerIndex{rhoR_s, rhoVR_s, BR_s, EtotR_s};
    //         auto const uL_ss = PerIndex{rhoL_s, rhoVL_ss, BL_ss, EtotL_ss};
    //         auto const uR_ss = PerIndex{rhoR_s, rhoVR_ss, BR_ss, EtotR_ss};
    //
    //         return std::make_tuple(uL_s, uL_ss, uR_ss, uR_s);
    //     };
    //
    //     return hlld_intermediate_states(fL, fR);
    // }

    template<auto direction>
    PHARE_ALWAYS_INLINE auto hlld_speeds_(auto const& uL, auto const& uR, auto const& jL,
                                          auto const& jR)
    {
        auto const BdotBL = uL.B.x * uL.B.x + uL.B.y * uL.B.y + uL.B.z * uL.B.z;
        auto const BdotBR = uR.B.x * uR.B.x + uR.B.y * uR.B.y + uR.B.z * uR.B.z;

        auto compute_hll_speeds
            = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR, auto VcompL,
                  auto VcompR, auto BcompL, auto BcompR) PHARE_LAMBDA_ALWAYS_INLINE {
                  auto cfastL = compute_fast_magnetosonic_(gamma_, uL.rho, BcompL, BdotBL, uL.P);
                  auto cfastR = compute_fast_magnetosonic_(gamma_, uR.rho, BcompR, BdotBR, uR.P);
                  auto SL     = std::min(VcompL - cfastL, VcompR - cfastR);
                  auto SR     = std::max(VcompL + cfastL, VcompR + cfastR);

                  return std::make_tuple(SL, SR);
              };

        auto compute_hlld_speeds
            = [&](auto rhoL, auto rhoR, auto PL, auto PR, auto BdotBL, auto BdotBR, auto VcompL,
                  auto VcompR, auto BcompL, auto BcompR) PHARE_LAMBDA_ALWAYS_INLINE {
                  auto [SL, SR] = compute_hll_speeds(rhoL, rhoR, PL, PR, BdotBL, BdotBR, VcompL,
                                                     VcompR, BcompL, BcompR);

                  auto PtL = PL + 0.5 * BdotBL;
                  auto PtR = PR + 0.5 * BdotBR;

                  auto const Bn = BcompL; // is the same on both sides

                  auto SM_numerator
                      = rhoR * VcompR * (SR - VcompR) - rhoL * VcompL * (SL - VcompL) - PtR + PtL;
                  auto SM_denominator = rhoR * (SR - VcompR) - rhoL * (SL - VcompL);
                  auto SM             = SM_numerator / SM_denominator;

                  auto const rhoL_s = rhoL * (SL - VcompL) / (SL - SM);
                  auto const rhoR_s = rhoR * (SR - VcompR) / (SR - SM);

                  auto const SL_s = SM - std::abs(Bn) / std::sqrt(rhoL_s);
                  auto const SR_s = SM + std::abs(Bn) / std::sqrt(rhoR_s);

                  auto hlld_speeds = std::make_tuple(SL, SL_s, SM, SR_s, SR);

                  uct_coefs_<direction>(uL, uR, jL, jR, hlld_speeds);


                  return hlld_speeds;
              };


        if constexpr (direction == Direction::X)
            return compute_hlld_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.x, uR.V.x,
                                       uL.B.x, uR.B.x);
        else if constexpr (direction == Direction::Y)
            return compute_hlld_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.y, uR.V.y,
                                       uL.B.y, uR.B.y);
        else if constexpr (direction == Direction::Z)
            return compute_hlld_speeds(uL.rho, uR.rho, uL.P, uR.P, BdotBL, BdotBR, uL.V.z, uR.V.z,
                                       uL.B.z, uR.B.z);
    }

    PHARE_ALWAYS_INLINE auto hlld_(auto const& uL, auto const& uL_s, auto const& uL_ss,
                                   auto const& uR_ss, auto const& uR_s, auto const& uR,
                                   auto const& fL, auto const& fR, auto const& hlld_speeds) const
    {
        auto const [SL, SL_s, SM, SR_s, SR] = hlld_speeds;

        auto const hlld = [&](auto const ul, auto const ul_s, auto const ul_ss, auto const ur_ss,
                              auto const ur_s, auto const ur, auto const fl,
                              auto const fr) PHARE_LAMBDA_ALWAYS_INLINE {
            auto const fl_s  = fl + SL * ul_s - SL * ul;
            auto const fr_s  = fr + SR * ur_s - SR * ur;
            auto const fl_ss = fl + SL_s * ul_ss - (SL_s - SL) * ul_s - SL * ul;
            auto const fr_ss = fr + SR_s * ur_ss - (SR_s - SR) * ur_s - SR * ur;

            // innermost region first, one two-way select at a time: a single six-way select is
            // one branch join that GCC cannot if-convert
            auto const f_ss = SM >= 0 ? fl_ss : fr_ss; // L** or R**
            auto const f_sr = SR_s <= 0 ? fr_s : f_ss; // R*
            auto const f_sl = SL_s >= 0 ? fl_s : f_sr; // L*
            auto const f_r  = SR <= 0.0 ? fr : f_sl;   // R
            return SL >= 0.0 ? fl : f_r;               // L
        };

        // one explicit call per component: a generic for_N here is left out of line in the large
        // face loop, which then stays scalar
        auto const flux = [&](auto const& get) PHARE_LAMBDA_ALWAYS_INLINE {
            return hlld(get(uL), get(uL_s), get(uL_ss), get(uR_ss), get(uR_s), get(uR), get(fL),
                        get(fR));
        };

        return PerIndex{flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.rho; }),
                        {flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.V.x; }),
                         flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.V.y; }),
                         flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.V.z; })},
                        {flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.B.x; }),
                         flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.B.y; }),
                         flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.B.z; })},
                        flux([](auto const& u) PHARE_LAMBDA_ALWAYS_INLINE { return u.P; })};
    }

    // in the hllc fallback used in idefix/pluto, they use hll averages for the magnetic field. we
    // do the same here for consistency.
    template<auto direction>
    PHARE_ALWAYS_INLINE void uct_coefs_(auto const& uL, auto const& uR, auto const& hlld_speeds)
    {
        auto const [SL, SL_s, SM, SR_s, SR] = hlld_speeds;

        // clipped speeds, as in hll.hpp
        SL_ = std::min(0.0, SL);
        SR_ = std::max(0.0, SR);

        bool const fallback
            = ((SL_s - SL) < 1.0e-4 * (SM - SL)) | ((SR_s - SR) > -1.0e-4 * (SR - SM));

        // both branches are computed and selected, so the face loop has no control flow
        auto const inv = 1.0 / (SR_ - SL_);

        auto const nuL = (SL_s + SL) / (std::abs(SL_s) + std::abs(SL));
        auto const nuR = (SR_s + SR) / (std::abs(SR_s) + std::abs(SR));

        auto const nu_s_ratio = (SR_s + SL_s) / (std::abs(SR_s) + std::abs(SL_s));
        auto const nu_s = std::abs(SR_s - SL_s) > 1.0e-9 * std::abs(SR - SL) ? nu_s_ratio : 0.0;

        auto const xiL = ((uL.V(direction) - SM) * (SL - SM)) / (SL_s + SL - 2. * SM);
        auto const xiR = ((uR.V(direction) - SM) * (SR - SM)) / (SR_s + SR - 2. * SM);

        auto const hll_d = -SR_ * SL_ * inv;

        uct_coefs[0] = fallback ? SR_ * inv : (1 + nu_s) * 0.5;
        uct_coefs[1] = fallback ? -SL_ * inv : (1 - nu_s) * 0.5;
        uct_coefs[2]
            = fallback ? hll_d : 0.5 * (nuL - nu_s) * xiL + 0.5 * (std::abs(SL_s) - nu_s * SL_s);
        uct_coefs[3]
            = fallback ? hll_d : 0.5 * (nuR - nu_s) * xiR + 0.5 * (std::abs(SR_s) - nu_s * SR_s);

        vt = vector_riemann_averaging(uL.V, uR.V);
    }

    template<auto direction>
    PHARE_ALWAYS_INLINE void uct_coefs_(auto const& uL, auto const& uR, auto const& jL,
                                        auto const& jR, auto const& hlld_speeds)
    {
        uct_coefs_<direction>(uL, uR, hlld_speeds);

        jt = vector_riemann_averaging(jL, jR);

        rhot = riemann_averaging(uL.rho, uR.rho);
    }
};
} // namespace PHARE::core

#endif
