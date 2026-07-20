#ifndef PHARE_CORE_NUMERICS_THERMO_THERMO_DEFS_HPP
#define PHARE_CORE_NUMERICS_THERMO_THERMO_DEFS_HPP

#include <stdexcept>
#include <string>
#include <unordered_map>

namespace PHARE::core
{
/**
 * @brief Enumerates all available thermodynamic (equation-of-state) models.
 *
 * Add a new enumerator here whenever a new concrete Thermo implementation is
 * introduced, and update getThermoModelFromString() accordingly.
 */
enum class ThermoModel {
    ideal_gas, ///< Ideal gas: P = ρ (γ - 1) u, constant heat capacity ratio γ.
};

/**
 * @brief Convert a string name to the corresponding ThermoModel enumerator.
 *
 * Expected string values (case-sensitive):
 *   - "ideal_gas"
 *
 * @param s  String representation of the thermodynamic model, typically read
 *           from the user configuration (e.g. \c dict["eos"].to<std::string>()).
 * @return   The matching ThermoModel enumerator.
 * @throws   std::runtime_error if @p s does not match any known model.
 */
inline ThermoModel getThermoModelFromString(std::string const& s)
{
    static std::unordered_map<std::string, ThermoModel> const map{
        {"ideal_gas", ThermoModel::ideal_gas},
    };

    auto const it = map.find(s);
    if (it == map.end())
        throw std::runtime_error("Unknown thermodynamic model: \"" + s + "\"");

    return it->second;
}

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_THERMO_THERMO_DEFS_HPP
