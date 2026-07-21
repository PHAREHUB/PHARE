#ifndef PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_CONVERSION_UTILS_HPP
#define PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_CONVERSION_UTILS_HPP

namespace PHARE::core
{
/**
 * @brief Compute the volumetric internal energy from the total energy.
 *
 * Strips the kinetic and magnetic contributions from the total energy:
 *
 *   e_int = Etot - ½ρ|v|² - ½|B|²
 *
 * This operation is EOS-agnostic: it depends only on the MHD energy
 * decomposition, not on any thermodynamic model.
 *
 * @param Etot  Total energy per unit volume.
 * @param rho   Mass density ρ.
 * @param vx    x-component of velocity.
 * @param vy    y-component of velocity.
 * @param vz    z-component of velocity.
 * @param bx    x-component of magnetic field.
 * @param by    y-component of magnetic field.
 * @param bz    z-component of magnetic field.
 * @return      Volumetric internal energy e_int = ρ u.
 */
inline auto internalEnergyFromTotalEnergy(auto Etot, auto rho,
                                          auto vx, auto vy, auto vz,
                                          auto bx, auto by, auto bz)
{
    return Etot - 0.5 * rho * (vx * vx + vy * vy + vz * vz)
                - 0.5 * (bx * bx + by * by + bz * bz);
}

/**
 * @brief Overload taking indexable velocity and magnetic field objects.
 *
 * @param Etot  Total energy per unit volume.
 * @param rho   Mass density ρ.
 * @param v     Indexable velocity, accessed as v[0], v[1], v[2].
 * @param B     Indexable magnetic field, accessed as B[0], B[1], B[2].
 * @return      Volumetric internal energy e_int = ρ u.
 */
inline auto internalEnergyFromTotalEnergy(auto Etot, auto rho, auto const& v, auto const& B)
{
    return internalEnergyFromTotalEnergy(Etot, rho, v[0], v[1], v[2], B[0], B[1], B[2]);
}

/**
 * @brief Compute the total energy from the volumetric internal energy.
 *
 * Adds the kinetic and magnetic contributions to the internal energy:
 *
 *   Etot = e_int + ½ρ|v|² + ½|B|²
 *
 * This operation is EOS-agnostic: it depends only on the MHD energy
 * decomposition, not on any thermodynamic model.  The caller is responsible
 * for obtaining @p e_int from the appropriate @c Thermo object.
 *
 * @param e_int Volumetric internal energy e_int = ρ u.
 * @param rho   Mass density ρ.
 * @param vx    x-component of velocity.
 * @param vy    y-component of velocity.
 * @param vz    z-component of velocity.
 * @param bx    x-component of magnetic field.
 * @param by    y-component of magnetic field.
 * @param bz    z-component of magnetic field.
 * @return      Total energy per unit volume Etot.
 */
inline auto totalEnergyFromInternalEnergy(auto e_int, auto rho,
                                          auto vx, auto vy, auto vz,
                                          auto bx, auto by, auto bz)
{
    return e_int + 0.5 * rho * (vx * vx + vy * vy + vz * vz)
                 + 0.5 * (bx * bx + by * by + bz * bz);
}

/**
 * @brief Overload taking indexable velocity and magnetic field objects.
 *
 * @param e_int Volumetric internal energy e_int = ρ u.
 * @param rho   Mass density ρ.
 * @param v     Indexable velocity, accessed as v[0], v[1], v[2].
 * @param B     Indexable magnetic field, accessed as B[0], B[1], B[2].
 * @return      Total energy per unit volume Etot.
 */
inline auto totalEnergyFromInternalEnergy(auto e_int, auto rho, auto const& v, auto const& B)
{
    return totalEnergyFromInternalEnergy(e_int, rho, v[0], v[1], v[2], B[0], B[1], B[2]);
}

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_CONVERSION_UTILS_HPP
