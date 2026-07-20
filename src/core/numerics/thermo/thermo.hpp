#ifndef PHARE_CORE_NUMERICS_THERMO_THERMO_HPP
#define PHARE_CORE_NUMERICS_THERMO_THERMO_HPP

#include "core/numerics/thermo/thermo_defs.hpp"

namespace PHARE::core
{
/**
 * @brief Abstract thermodynamic equation-of-state interface.
 *
 * Inspired by Cantera's ThermoPhase, this class defines the contract that any
 * equation of state (EOS) must satisfy. The thermodynamic state is set through
 * one of the provided setters, after which all derived quantities can be queried.
 *
 * The object holds the state internally; callers set it once and query as needed,
 * which amortises any conversion cost when multiple quantities are required.
 *
 * @note Total energy in MHD includes kinetic (½ρv²) and magnetic (½B²) contributions
 *       that are EOS-agnostic. This interface covers only the thermodynamic (internal
 *       energy) part; callers are responsible for adding those external contributions.
 *
 * @note The EOS type is a runtime user choice (read from the configuration dictionary).
 *       Concrete implementations are constructed via a factory and owned as
 *       std::unique_ptr<Thermo> by the solver or state object.
 */
class Thermo
{
public:
    virtual ~Thermo() = default;

    /**
     * @brief Set the thermodynamic state from mass density and pressure.
     *
     * This is the natural setter for the primitive-variable path, where ρ and P
     * are directly available. The specific internal energy is derived from these
     * two quantities via the EOS.
     *
     * @param rho Mass density ρ (kg/m³ in SI, or simulation units).
     * @param p   Thermal pressure P.
     */
    virtual void setState_DP(double rho, double p) = 0;

    /**
     * @brief Set the thermodynamic state from specific internal energy and pressure.
     *
     * This is the natural setter for the conservative-variable path: after stripping
     * kinetic and magnetic contributions from total energy, the caller is left with
     * the volumetric internal energy e = Etot - ½ρv² - ½B², from which the specific
     * internal energy is u = e/ρ. Together with pressure, these two quantities uniquely
     * determine the full thermodynamic state (e.g. for an ideal gas, density is recovered
     * as ρ = P / ((γ-1)·u)).
     *
     * @param u   Specific internal energy per unit mass u (J/kg in SI, or simulation units).
     * @param p   Thermal pressure P.
     */
    virtual void setState_UP(double u, double p) = 0;

    /**
     * @brief Set the thermodynamic state from mass density and specific internal energy.
     *
     * This is the natural setter when recovering pressure from the conservative variables:
     * after stripping kinetic and magnetic contributions from total energy, the caller has
     * the volumetric internal energy e = Etot - ½ρv² - ½B² and the density ρ, which
     * together fully determine the thermodynamic state via the EOS (e.g. for an ideal gas,
     * P = ρ (γ-1) u).
     *
     * @param rho Mass density ρ.
     * @param u   Specific internal energy per unit mass u (J/kg in SI, or simulation units).
     */
    virtual void setState_DU(double rho, double u) = 0;

    /**
     * @brief Return the thermal pressure.
     * @pre State must have been set via one of the setState_* methods.
     * @return Thermal pressure P.
     */
    virtual double pressure() const = 0;

    /**
     * @brief Return the specific internal energy (internal energy per unit mass).
     * @pre State must have been set via one of the setState_* methods.
     * @return Specific internal energy u.
     */
    virtual double internalEnergy() const = 0;

    /**
     * @brief Return the adiabatic sound speed.
     * @pre State must have been set via one of the setState_* methods.
     * @return Sound speed cs = √(∂P/∂ρ|_s).
     */
    virtual double soundSpeed() const = 0;

    /**
     * @brief Return the temperature.
     * @pre State must have been set via one of the setState_* methods.
     * @return Temperature T (in simulation units; for ideal gas T = P/ρ in normalised units).
     */
    virtual double temperature() const = 0;

    /**
     * @brief Return the thermodynamic model implemented by this object.
     * @return A ThermoModel enumerator identifying the concrete EOS.
     */
    virtual ThermoModel type() const = 0;
};
} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_THERMO_THERMO_HPP
