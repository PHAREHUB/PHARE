#ifndef CORE_NUMERICS_MHD_EQUATIONS_HPP
#define CORE_NUMERICS_MHD_EQUATIONS_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"

// the magnetic fluxes computations should be removed from here
namespace PHARE::core
{
template<bool Hall>
class MHDEquations
{
public:
    constexpr static bool hall = Hall;

    MHDEquations(double const gamma)
        : gamma_{gamma}
    {
    }

    template<auto direction>
    auto compute(auto const& u) const
    {
        auto const rho = u.rho;
        auto const V   = u.V;
        auto const B   = u.totalB();
        auto const B1  = u.B1;
        auto const B0  = u.B0;
        auto const P   = u.P;

        auto const MagPressure = 0.5 * (B1.x * B1.x + B1.y * B1.y + B1.z * B1.z)
                                 + (B1.x * B0.x + B1.y * B0.y + B1.z * B0.z);
        auto const GeneralisedPressure = P + MagPressure;

        auto const HydroEnergy = P / (gamma_ - 1.0) + kineticEnergy(rho, V.x, V.y, V.z);
        auto const Ex          = V.z * B.y - V.y * B.z;
        auto const Ey          = V.x * B.z - V.z * B.x;
        auto const Ez          = V.y * B.x - V.x * B.y;

        if constexpr (direction == Direction::X)
        {
            auto F_rho   = rho * V.x;
            auto F_rhoVx
                = rho * V.x * V.x + GeneralisedPressure - (B1.x * B1.x + B1.x * B0.x + B0.x * B1.x);
            auto F_rhoVy = rho * V.x * V.y - (B1.x * B1.y + B1.x * B0.y + B0.x * B1.y);
            auto F_rhoVz = rho * V.x * V.z - (B1.x * B1.z + B1.x * B0.z + B0.x * B1.z);
            auto F_Bx    = 0.0;
            auto F_By    = B.y * V.x - V.y * B.x;
            auto F_Bz    = B.z * V.x - V.z * B.x;
            auto F_Etot  = (HydroEnergy + P) * V.x + Ey * B1.z - Ez * B1.y;

            return PerIndex{F_rho, {F_rhoVx, F_rhoVy, F_rhoVz}, {F_Bx, F_By, F_Bz}, F_Etot};
        }
        if constexpr (direction == Direction::Y)
        {
            auto F_rho   = rho * V.y;
            auto F_rhoVx = rho * V.y * V.x - (B1.y * B1.x + B1.y * B0.x + B0.y * B1.x);
            auto F_rhoVy
                = rho * V.y * V.y + GeneralisedPressure - (B1.y * B1.y + B1.y * B0.y + B0.y * B1.y);
            auto F_rhoVz = rho * V.y * V.z - (B1.y * B1.z + B1.y * B0.z + B0.y * B1.z);
            auto F_Bx    = B.x * V.y - V.x * B.y;
            auto F_By    = 0.0;
            auto F_Bz    = B.z * V.y - V.z * B.y;
            auto F_Etot  = (HydroEnergy + P) * V.y + Ez * B1.x - Ex * B1.z;

            return PerIndex{F_rho, {F_rhoVx, F_rhoVy, F_rhoVz}, {F_Bx, F_By, F_Bz}, F_Etot};
        }
        if constexpr (direction == Direction::Z)
        {
            auto F_rho   = rho * V.z;
            auto F_rhoVx = rho * V.z * V.x - (B1.z * B1.x + B1.z * B0.x + B0.z * B1.x);
            auto F_rhoVy = rho * V.z * V.y - (B1.z * B1.y + B1.z * B0.y + B0.z * B1.y);
            auto F_rhoVz
                = rho * V.z * V.z + GeneralisedPressure - (B1.z * B1.z + B1.z * B0.z + B0.z * B1.z);
            auto F_Bx    = B.x * V.z - V.x * B.z;
            auto F_By    = B.y * V.z - V.y * B.z;
            auto F_Bz    = 0.0;
            auto F_Etot  = (HydroEnergy + P) * V.z + Ex * B1.y - Ey * B1.x;

            return PerIndex{F_rho, {F_rhoVx, F_rhoVy, F_rhoVz}, {F_Bx, F_By, F_Bz}, F_Etot};
        }
    }

    template<auto direction>
    auto compute(auto const& u, auto const& J) const
    {
        PerIndex f = compute<direction>(u);

        if constexpr (Hall)
            hall_contribution_<direction>(u.rho, u.totalB(), u.B1, J, f.B1, f.P);

        return f;
    }

    template<auto direction>
    void dissipative_contributions(auto const& Et, auto const& Bt, auto& F_B, auto& F_Etot) const
    {
        if constexpr (direction == Direction::X)
        {
            F_B.y += -Et.z;
            F_B.z += Et.y;
            F_Etot += Et.y * Bt.z - Et.z * Bt.y;
        }
        if constexpr (direction == Direction::Y)
        {
            F_B.x += Et.z;
            F_B.z += -Et.x;
            F_Etot += Et.z * Bt.x - Et.x * Bt.z;
        }
        if constexpr (direction == Direction::Z)
        {
            F_B.x += -Et.y;
            F_B.y += Et.x;
            F_Etot += Et.x * Bt.y - Et.y * Bt.x;
        }
    }

private:
    double const gamma_;

    template<auto direction>
    void hall_contribution_(auto const& rho, auto const& B, auto const& B1, auto const& J,
                            auto& F_B, auto& F_Etot) const
    {
        auto const invRho = 1.0 / rho;

        auto const JxB_x = J.y * B.z - J.z * B.y;
        auto const JxB_y = J.z * B.x - J.x * B.z;
        auto const JxB_z = J.x * B.y - J.y * B.x;

        auto const B1dotJ = B1.x * J.x + B1.y * J.y + B1.z * J.z;
        auto const BdotB1 = B.x * B1.x + B.y * B1.y + B.z * B1.z;

        if constexpr (direction == Direction::X)
        {
            F_B.y += -JxB_z * invRho;
            F_B.z += JxB_y * invRho;
            F_Etot += (B1dotJ * B.x - BdotB1 * J.x) * invRho;
        }
        if constexpr (direction == Direction::Y)
        {
            F_B.x += JxB_z * invRho;
            F_B.z += -JxB_x * invRho;
            F_Etot += (B1dotJ * B.y - BdotB1 * J.y) * invRho;
        }
        if constexpr (direction == Direction::Z)
        {
            F_B.x += -JxB_y * invRho;
            F_B.y += JxB_x * invRho;
            F_Etot += (B1dotJ * B.z - BdotB1 * J.z) * invRho;
        }
    }
};

} // namespace PHARE::core

#endif
