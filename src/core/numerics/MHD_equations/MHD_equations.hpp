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
    auto compute(auto const& u) const
    {
        auto const rho = u.rho;
        auto const V   = u.V;
        auto const B   = u.B;
        auto const P   = u.P;

        auto const GeneralisedPressure = P + 0.5 * (B.x * B.x + B.y * B.y + B.z * B.z);
        auto const TotalEnergy         = eosPToEtot(gamma_, rho, V.x, V.y, V.z, B.x, B.y, B.z, P);

        if constexpr (direction == Direction::X)
        {
            auto F_rho   = rho * V.x;
            auto F_rhoVx = rho * V.x * V.x + GeneralisedPressure - B.x * B.x;
            auto F_rhoVy = rho * V.x * V.y - B.x * B.y;
            auto F_rhoVz = rho * V.x * V.z - B.x * B.z;
            auto F_Bx    = 0.0;
            auto F_By    = B.y * V.x - V.y * B.x;
            auto F_Bz    = B.z * V.x - V.z * B.x;
            auto F_Etot  = (TotalEnergy + GeneralisedPressure) * V.x
                          - B.x * (V.x * B.x + V.y * B.y + V.z * B.z);

            return PerIndex{F_rho, {F_rhoVx, F_rhoVy, F_rhoVz}, {F_Bx, F_By, F_Bz}, F_Etot};
        }
        if constexpr (direction == Direction::Y)
        {
            auto F_rho   = rho * V.y;
            auto F_rhoVx = rho * V.y * V.x - B.y * B.x;
            auto F_rhoVy = rho * V.y * V.y + GeneralisedPressure - B.y * B.y;
            auto F_rhoVz = rho * V.y * V.z - B.y * B.z;
            auto F_Bx    = B.x * V.y - V.x * B.y;
            auto F_By    = 0.0;
            auto F_Bz    = B.z * V.y - V.z * B.y;
            auto F_Etot  = (TotalEnergy + GeneralisedPressure) * V.y
                          - B.y * (V.x * B.x + V.y * B.y + V.z * B.z);

            return PerIndex{F_rho, {F_rhoVx, F_rhoVy, F_rhoVz}, {F_Bx, F_By, F_Bz}, F_Etot};
        }
        if constexpr (direction == Direction::Z)
        {
            auto F_rho   = rho * V.z;
            auto F_rhoVx = rho * V.z * V.x - B.z * B.x;
            auto F_rhoVy = rho * V.z * V.y - B.z * B.y;
            auto F_rhoVz = rho * V.z * V.z + GeneralisedPressure - B.z * B.z;
            auto F_Bx    = B.x * V.z - V.x * B.z;
            auto F_By    = B.y * V.z - V.y * B.z;
            auto F_Bz    = 0.0;
            auto F_Etot  = (TotalEnergy + GeneralisedPressure) * V.z
                          - B.z * (V.x * B.x + V.y * B.y + V.z * B.z);

            return PerIndex{F_rho, {F_rhoVx, F_rhoVy, F_rhoVz}, {F_Bx, F_By, F_Bz}, F_Etot};
        }
    }

    template<auto direction>
    auto compute(auto const& u, auto const& J) const
    {
        PerIndex f = compute<direction>(u);

        if constexpr (Hall)
            hall_contribution_<direction>(u.rho, u.B, J, f.B, f.P);
        if constexpr (Resistivity)
            resistive_contributions_<direction>(eta_, u.B, J, f.B, f.P);

        return f;
    }

    template<auto direction>
    auto compute(auto const& u, auto const& J, auto const& LaplJ) const
    {
        PerIndex f = compute<direction>(u);

        if constexpr (Hall)
            hall_contribution_<direction>(u.rho, u.B, J, f.B, f.P);
        if constexpr (Resistivity)
            resistive_contributions_<direction>(eta_, u.B, J, f.B, f.P);

        resistive_contributions_<direction>(nu_, u.B, J, f.B, f.P);

        return f;
    }


private:
    double const gamma_;
    double const eta_;
    double const nu_;

    template<auto direction>
    void hall_contribution_(auto const& rho, auto const& B, auto const& J, auto& F_B,
                            auto& F_Etot) const
    {
        auto const invRho = 1.0 / rho;

        auto const JxB_x = J.y * B.z - J.z * B.y;
        auto const JxB_y = J.z * B.x - J.x * B.z;
        auto const JxB_z = J.x * B.y - J.y * B.x;

        auto const BdotJ = B.x * J.x + B.y * J.y + B.z * J.z;
        auto const BdotB = B.x * B.x + B.y * B.y + B.z * B.z;

        if constexpr (direction == Direction::X)
        {
            F_B.y += -JxB_z * invRho;
            F_B.z += JxB_y * invRho;
            F_Etot += (BdotJ * B.x - BdotB * J.x) * invRho;
        }
        if constexpr (direction == Direction::Y)
        {
            F_B.x += JxB_z * invRho;
            F_B.z += -JxB_x * invRho;
            F_Etot += (BdotJ * B.y - BdotB * J.y) * invRho;
        }
        if constexpr (direction == Direction::Z)
        {
            F_B.x += -JxB_y * invRho;
            F_B.y += JxB_x * invRho;
            F_Etot += (BdotJ * B.z - BdotB * J.z) * invRho;
        }
    }

    template<auto direction>
    void resistive_contributions_(auto const& coef, auto const& B, auto const& J, auto& F_B,
                                  auto& F_Etot) const
    // Can be used for both resistivity with J and eta and hyper resistivity with laplJ and nu
    {
        if constexpr (direction == Direction::X)
        {
            F_B.y += -J.z * coef;
            F_B.z += J.y * coef;
            F_Etot += (J.y * B.z - J.z * B.y) * coef;
        }
        if constexpr (direction == Direction::Y)
        {
            F_B.x += J.z * coef;
            F_B.z += -J.x * coef;
            F_Etot += (J.z * B.x - J.x * B.z) * coef;
        }
        if constexpr (direction == Direction::Z)
        {
            F_B.x += -J.y * coef;
            F_B.y += J.x * coef;
            F_Etot += (J.x * B.y - J.y * B.x) * coef;
        }
    }
};

} // namespace PHARE::core

#endif
