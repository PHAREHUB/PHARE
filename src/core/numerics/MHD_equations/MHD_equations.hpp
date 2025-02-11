#ifndef CORE_NUMERICS_MHD_EQUATIONS_HPP
#define CORE_NUMERICS_MHD_EQUATIONS_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"

namespace PHARE::core
{
template<bool Hall, bool Resistivity, bool HyperResistivity>
class MHDEquations
{
public:
    constexpr static bool hall             = Hall;
    constexpr static bool resistivity      = Resistivity;
    constexpr static bool hyperResistivity = HyperResistivity;

    MHDEquations(double const gamma, double const eta, double const nu)
        : gamma_{gamma}
        , eta_{eta}
        , nu_{nu}
    {
    }

    template<auto direction>
    PerIndex compute(PerIndex const& u) const
    {
        auto const rho = u.rho;
        auto const Vx  = u.Vx;
        auto const Vy  = u.Vy;
        auto const Vz  = u.Vz;
        auto const Bx  = u.Bx;
        auto const By  = u.By;
        auto const Bz  = u.Bz;
        auto const P   = u.P;

        auto const GeneralisedPressure = P + 0.5 * (Bx * Bx + By * By + Bz * Bz);
        auto const TotalEnergy         = eosPToEtot(gamma_, rho, Vx, Vy, Vz, Bx, By, Bz, P);

        if constexpr (direction == Direction::X)
        {
            auto F_rho   = rho * Vx;
            auto F_rhoVx = rho * Vx * Vx + GeneralisedPressure - Bx * Bx;
            auto F_rhoVy = rho * Vx * Vy - Bx * By;
            auto F_rhoVz = rho * Vx * Vz - Bx * Bz;
            auto F_Bx    = 0.0;
            auto F_By    = By * Vx - Vy * Bx;
            auto F_Bz    = Bz * Vx - Vz * Bx;
            auto F_Etot
                = (TotalEnergy + GeneralisedPressure) * Vx - Bx * (Vx * Bx + Vy * By + Vz * Bz);

            return PerIndex(F_rho, F_rhoVx, F_rhoVy, F_rhoVz, F_Bx, F_By, F_Bz, F_Etot);
        }
        if constexpr (direction == Direction::Y)
        {
            auto F_rho   = rho * Vy;
            auto F_rhoVx = rho * Vy * Vx - By * Bx;
            auto F_rhoVy = rho * Vy * Vy + GeneralisedPressure - By * By;
            auto F_rhoVz = rho * Vy * Vz - By * Bz;
            auto F_Bx    = Bx * Vy - Vx * By;
            auto F_By    = 0.0;
            auto F_Bz    = Bz * Vy - Vz * By;
            auto F_Etot
                = (TotalEnergy + GeneralisedPressure) * Vy - By * (Vx * Bx + Vy * By + Vz * Bz);

            return PerIndex(F_rho, F_rhoVx, F_rhoVy, F_rhoVz, F_Bx, F_By, F_Bz, F_Etot);
        }
        if constexpr (direction == Direction::Z)
        {
            auto F_rho   = rho * Vz;
            auto F_rhoVx = rho * Vz * Vx - Bz * Bx;
            auto F_rhoVy = rho * Vz * Vy - Bz * By;
            auto F_rhoVz = rho * Vz * Vz + GeneralisedPressure - Bz * Bz;
            auto F_Bx    = Bx * Vz - Vx * Bz;
            auto F_By    = By * Vz - Vy * Bz;
            auto F_Bz    = 0.0;
            auto F_Etot
                = (TotalEnergy + GeneralisedPressure) * Vz - Bz * (Vx * Bx + Vy * By + Vz * Bz);

            return PerIndex(F_rho, F_rhoVx, F_rhoVy, F_rhoVz, F_Bx, F_By, F_Bz, F_Etot);
        }
    }

    template<auto direction>
    PerIndex compute(PerIndex const& u, auto const& j) const
    {
        PerIndex f               = compute<direction>(u);
        auto const& [jx, jy, jz] = j;

        if constexpr (Hall)
            hall_contribution_<direction>(u.rho, u.Bx, u.By, u.Bz, jx, jy, jz, f.Bx, f.By, f.Bz,
                                          f.P);
        if constexpr (Resistivity)
            resistive_contributions_<direction>(eta_, u.Bx, u.By, u.Bz, jx, jy, jz, f.Bx, f.By,
                                                f.Bz, f.P);

        return f;
    }

    template<auto direction>
    PerIndex compute(PerIndex const& u, auto const& j, auto const& LaplJ) const
    {
        PerIndex f               = compute<direction>(u);
        auto const& [jx, jy, jz] = j;

        if constexpr (Hall)
            hall_contribution_<direction>(u.rho, u.Bx, u.By, u.Bz, jx, jy, jz, f.Bx, f.By, f.Bz,
                                          f.P);
        if constexpr (Resistivity)
            resistive_contributions_<direction>(eta_, u.Bx, u.By, u.Bz, jx, jy, jz, f.Bx, f.By,
                                                f.Bz, f.P);

        resistive_contributions_<direction>(nu_, u.Bx, u.By, u.Bz, jx, jy, jz, f.Bx, f.By, f.Bz,
                                            f.P);

        return f;
    }


private:
    double const gamma_;
    double const eta_;
    double const nu_;

    template<auto direction>
    void hall_contribution_(auto const& rho, auto const& Bx, auto const& By, auto const& Bz,
                            auto const& Jx, auto const& Jy, auto const& Jz, auto& F_Bx, auto& F_By,
                            auto& F_Bz, auto& F_Etot) const
    {
        auto invRho = 1.0 / rho;

        auto JxB_x = Jy * Bz - Jz * By;
        auto JxB_y = Jz * Bx - Jx * Bz;
        auto JxB_z = Jx * By - Jy * Bx;

        auto BdotJ = Bx * Jx + By * Jy + Bz * Jz;
        auto BdotB = Bx * Bx + By * By + Bz * Bz;

        if constexpr (direction == Direction::X)
        {
            F_By += -JxB_z * invRho;
            F_Bz += JxB_y * invRho;
            F_Etot += (BdotJ * Bx - BdotB * Jx) * invRho;
        }
        if constexpr (direction == Direction::Y)
        {
            F_Bx += JxB_z * invRho;
            F_Bz += -JxB_x * invRho;
            F_Etot += (BdotJ * By - BdotB * Jy) * invRho;
        }
        if constexpr (direction == Direction::Z)
        {
            F_Bx += -JxB_y * invRho;
            F_By += JxB_x * invRho;
            F_Etot += (BdotJ * Bz - BdotB * Jz) * invRho;
        }
    }

    template<auto direction>
    void resistive_contributions_(auto const& pc, auto const& Bx, auto const& By, auto const& Bz,
                                  auto const& Jx, auto const& Jy, auto const& Jz, auto& F_Bx,
                                  auto& F_By, auto& F_Bz, auto& F_Etot) const
    // Can be used for both resistivity with J and eta and hyper resistivity with laplJ and nu
    {
        if constexpr (direction == Direction::X)
        {
            F_By += -Jz * pc;
            F_Bz += Jy * pc;
            F_Etot += (Jy * Bz - Jz * By) * pc;
        }
        if constexpr (direction == Direction::Y)
        {
            F_Bx += Jz * pc;
            F_Bz += -Jx * pc;
            F_Etot += (Jz * Bx - Jx * Bz) * pc;
        }
        if constexpr (direction == Direction::Y)
        {
            F_Bx += -Jy * pc;
            F_By += Jx * pc;
            F_Etot += (Jx * By - Jy * Bx) * pc;
        }
    }
};

} // namespace PHARE::core

#endif
