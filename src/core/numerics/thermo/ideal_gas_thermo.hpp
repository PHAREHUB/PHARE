#ifndef PHARE_CORE_NUMERICS_THERMO_IDEAL_GAS_THERMO_HPP
#define PHARE_CORE_NUMERICS_THERMO_IDEAL_GAS_THERMO_HPP

#include <cmath>

#include "core/numerics/thermo/thermo.hpp"
#include "core/numerics/thermo/thermo_defs.hpp"

namespace PHARE::core
{
/**
 * @brief Ideal gas equation of state.
 *
 * Implements the Thermo interface for a perfect gas with a constant heat
 * capacity ratio γ. The caloric equation of state is:
 *
 *   P = ρ (γ - 1) u
 *
 * where u is the specific internal energy (per unit mass).
 *
 * The canonical internal state is (ρ, u). Both setters convert their inputs
 * to this pair; all getters derive their return values from it.
 *
 * @note In the normalized units used by PHARE, temperature
 *       reduces to T = P / ρ = (γ - 1) u.
 */
class IdealGasThermo final : public Thermo
{
public:
    /**
     * @brief Construct an ideal gas EOS with the given heat capacity ratio.
     * @param gamma Heat capacity ratio γ = C_p / C_v  (γ > 1).
     */
    explicit IdealGasThermo(double gamma)
        : gamma_{gamma}
    {
    }

    /**
     * @brief Set the state from mass density and pressure.
     *
     * Derives the specific internal energy as u = P / (ρ (γ - 1)).
     *
     * @param rho Mass density ρ.
     * @param p   Thermal pressure P.
     */
    void setState_DP(double rho, double p) override
    {
        rho_ = rho;
        u_   = p / (rho * (gamma_ - 1.0));
    }

    /**
     * @brief Set the state from specific internal energy and pressure.
     *
     * Derives the mass density as ρ = P / ((γ - 1) u).
     *
     * @param u  Specific internal energy per unit mass u.
     * @param p  Thermal pressure P.
     */
    void setState_UP(double u, double p) override
    {
        u_   = u;
        rho_ = p / ((gamma_ - 1.0) * u);
    }

    /**
     * @brief Set the state from mass density and specific internal energy.
     *
     * Stores ρ and u directly; pressure is recovered as P = ρ (γ - 1) u.
     *
     * @param rho Mass density ρ.
     * @param u   Specific internal energy per unit mass u.
     */
    void setState_DU(double rho, double u) override
    {
        rho_ = rho;
        u_   = u;
    }

    /**
     * @brief Return the thermal pressure.
     * @return P = ρ (γ - 1) u.
     */
    double pressure() const override { return rho_ * (gamma_ - 1.0) * u_; }

    /**
     * @brief Return the specific internal energy (per unit mass).
     * @return u.
     */
    double internalEnergy() const override { return u_; }

    /**
     * @brief Return the adiabatic sound speed.
     * @return cs = √(γ P / ρ) = √(γ (γ - 1) u).
     */
    double soundSpeed() const override { return std::sqrt(gamma_ * pressure() / rho_); }

    /**
     * @brief Return the temperature in normalized units (k_B / m = 1).
     * @return T = P / ρ = (γ - 1) u.
     */
    double temperature() const override { return pressure() / rho_; }

    /**
     * @brief Return the thermodynamic model identifier.
     * @return ThermoModel::ideal_gas.
     */
    ThermoModel type() const override { return ThermoModel::ideal_gas; }

private:
    double gamma_;
    double rho_ = 0.;
    double u_   = 0.;
};

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_THERMO_IDEAL_GAS_THERMO_HPP
