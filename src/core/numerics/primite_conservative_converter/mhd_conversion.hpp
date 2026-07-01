#ifndef PHARE_CORE_NUMERICS_MHD_CONVERSION_HPP
#define PHARE_CORE_NUMERICS_MHD_CONVERSION_HPP

#include <tuple>

namespace PHARE::core
{

inline auto magneticEnergy(auto const& bx, auto const& by, auto const& bz)
{
    return 0.5 * (bx * bx + by * by + bz * bz);
}

inline auto kineticEnergy(auto const& rho, auto const& vx, auto const& vy, auto const& vz)
{
    return 0.5 * rho * (vx * vx + vy * vy + vz * vz);
}

inline auto totalMagneticComponents(auto const& b1x, auto const& b1y, auto const& b1z,
                                    auto const& b0x, auto const& b0y, auto const& b0z)
{
    return std::make_tuple(b1x + b0x, b1y + b0y, b1z + b0z);
}

inline auto totalMagneticEnergy(auto const& b1x, auto const& b1y, auto const& b1z, auto const& b0x,
                                auto const& b0y, auto const& b0z)
{
    auto const& [bx, by, bz] = totalMagneticComponents(b1x, b1y, b1z, b0x, b0y, b0z);
    return magneticEnergy(bx, by, bz);
}

// Etot1 (B1-only energy) -> Etot (total-field energy): add the B0 self-energy and the B0.B1
// cross term so that 1/2|B0+B1|^2 = 1/2|B1|^2 + B0.B1 + 1/2|B0|^2 is recovered. Used for
// diagnostics, which expose the total energy.
inline auto etot1ToEtot(auto const& etot1, auto const& b1x, auto const& b1y, auto const& b1z,
                        auto const& b0x, auto const& b0y, auto const& b0z)
{
    return etot1 + magneticEnergy(b0x, b0y, b0z) + (b0x * b1x + b0y * b1y + b0z * b1z);
}

inline auto etotToEtot1(auto const& etot, auto const& b1x, auto const& b1y, auto const& b1z,
                        auto const& b0x, auto const& b0y, auto const& b0z)
{
    return etot - magneticEnergy(b0x, b0y, b0z) - (b0x * b1x + b0y * b1y + b0z * b1z);
}

// Total-field equations of state (kept for diagnostics / classical use).
inline auto eosPToEtot(double const gamma, auto const& rho, auto const& vx, auto const& vy,
                       auto const& vz, auto const& bx, auto const& by, auto const& bz,
                       auto const& p)
{
    return p / (gamma - 1.0) + kineticEnergy(rho, vx, vy, vz) + magneticEnergy(bx, by, bz);
}

inline auto eosEtotToP(double const gamma, auto const& rho, auto const& vx, auto const& vy,
                       auto const& vz, auto const& bx, auto const& by, auto const& bz,
                       auto const& etot)
{
    return (gamma - 1.0) * (etot - kineticEnergy(rho, vx, vy, vz) - magneticEnergy(bx, by, bz));
}

// Split-field equations of state: the conserved energy Etot1 stores only the B1 magnetic
// energy, so the pressure solve uses B1 alone (B0 is a static background, not part of the
// evolving pressure).
inline auto eosPToEtot1(double const gamma, auto const& rho, auto const& vx, auto const& vy,
                        auto const& vz, auto const& b1x, auto const& b1y, auto const& b1z,
                        auto const& p)
{
    return p / (gamma - 1.0) + kineticEnergy(rho, vx, vy, vz) + magneticEnergy(b1x, b1y, b1z);
}

inline auto eosEtot1ToP(double const gamma, auto const& rho, auto const& vx, auto const& vy,
                        auto const& vz, auto const& b1x, auto const& b1y, auto const& b1z,
                        auto const& etot1)
{
    return (gamma - 1.0) * (etot1 - kineticEnergy(rho, vx, vy, vz) - magneticEnergy(b1x, b1y, b1z));
}

} // namespace PHARE::core

#endif
