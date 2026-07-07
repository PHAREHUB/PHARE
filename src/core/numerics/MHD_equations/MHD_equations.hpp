#ifndef CORE_NUMERICS_MHD_EQUATIONS_HPP
#define CORE_NUMERICS_MHD_EQUATIONS_HPP

#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"

// the magnetic fluxes computations should be removed from here
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
        auto const B   = u.totalB(); // total field, used by the induction and energy fluxes
        auto const B1  = u.B1;       // perturbation field (the evolved/reconstructed field)
        auto const B0  = u.B0;       // static background field
        auto const P   = u.P;

        // Well-balanced B-split Maxwell stress: total-field stress minus the static B0 self-stress
        // leaves only B1 and B1.B0 cross terms. Makes a B1=0/V=0/uniform-P state a steady state.
        auto const MagPressure = 0.5 * (B1.x * B1.x + B1.y * B1.y + B1.z * B1.z)
                                 + (B1.x * B0.x + B1.y * B0.y + B1.z * B0.z);
        auto const GeneralisedPressure = P + MagPressure;

        // Hydrodynamic energy (kinetic + thermal). Magnetic energy transport is the Poynting term
        // E x B1 below, with the motional field E = V x B (total B).
        auto const HydroEnergy = P / (gamma_ - 1.0) + 0.5 * rho * (V.x * V.x + V.y * V.y + V.z * V.z);
        auto const Ex          = V.z * B.y - V.y * B.z;
        auto const Ey          = V.x * B.z - V.z * B.x;
        auto const Ez          = V.y * B.x - V.x * B.y;

        if constexpr (direction == Direction::X)
        {
            auto F_rho   = rho * V.x;
            auto F_rhoVx = rho * V.x * V.x + GeneralisedPressure
                           - (B1.x * B1.x + B1.x * B0.x + B0.x * B1.x);
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
            auto F_rhoVy = rho * V.y * V.y + GeneralisedPressure
                           - (B1.y * B1.y + B1.y * B0.y + B0.y * B1.y);
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
            auto F_rhoVz = rho * V.z * V.z + GeneralisedPressure
                           - (B1.z * B1.z + B1.z * B0.z + B0.z * B1.z);
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
            hall_contribution_<direction>(u.rho, u.totalB(), J, f.B1, f.P);
        // if constexpr (Resistivity)
        //     resistive_contributions_<direction>(eta_, u.totalB(), J, f.B1, f.P);

        return f;
    }

    // template<auto direction>
    // auto compute(auto const& u, auto const& J, auto const& LaplJ) const
    // {
    //     PerIndex f = compute<direction>(u);
    //
    //     if constexpr (Hall)
    //         hall_contribution_<direction>(u.rho, u.B, J, f.B, f.P);
    //     if constexpr (Resistivity)
    //         resistive_contributions_<direction>(eta_, u.B, J, f.B, f.P);
    //
    //     resistive_contributions_<direction>(nu_, u.B, -LaplJ, f.B, f.P);
    //
    //     return f;
    // }

    template<auto direction>
    void resistive_contributions(auto const& coef, auto const& Bt, auto const& Jt, auto& F_B,
                                 auto& F_Etot) const
    // Can be used for both resistivity with J and eta and hyper resistivity with laplJ and nu. The
    // work is done on the tranverse riemann averaged components avoid extra reconstructions. This
    // optimisation is possible since these operations are linear.
    {
        if constexpr (direction == Direction::X)
        {
            F_B.y += -Jt.z * coef;
            F_B.z += Jt.y * coef;
            F_Etot += (Jt.y * Bt.z - Jt.z * Bt.y) * coef;
        }
        if constexpr (direction == Direction::Y)
        {
            F_B.x += Jt.z * coef;
            F_B.z += -Jt.x * coef;
            F_Etot += (Jt.z * Bt.x - Jt.x * Bt.z) * coef;
        }
        if constexpr (direction == Direction::Z)
        {
            F_B.x += -Jt.y * coef;
            F_B.y += Jt.x * coef;
            F_Etot += (Jt.x * Bt.y - Jt.y * Bt.x) * coef;
        }
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
};

} // namespace PHARE::core

#endif
