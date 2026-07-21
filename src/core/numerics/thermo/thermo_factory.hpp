#ifndef PHARE_CORE_NUMERICS_THERMO_THERMO_FACTORY_HPP
#define PHARE_CORE_NUMERICS_THERMO_THERMO_FACTORY_HPP

#include <memory>
#include <stdexcept>

#include "core/numerics/thermo/thermo.hpp"
#include "core/numerics/thermo/thermo_defs.hpp"
#include "core/numerics/thermo/ideal_gas_thermo.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::core
{
/**
 * @brief Construct a Thermo instance from the MHD state configuration dict.
 *
 * Reads @c dict["eos"] to select the concrete EOS model, then reads any
 * model-specific parameters (e.g. @c dict["heat_capacity_ratio"] for the
 * ideal gas). The returned object is heap-allocated and owned by the caller
 * via a shared pointer.
 *
 * Expected dict keys (under the @c mhd_state sub-dict):
 *   - @c "eos"                  : string, e.g. @c "ideal_gas"
 *   - @c "heat_capacity_ratio"  : double, required when @c eos == @c "ideal_gas"
 *
 * @param dict  The @c mhd_state sub-dict from the simulation PHAREDict.
 * @return      A shared pointer to the constructed Thermo object.
 * @throws      std::runtime_error if the EOS string is unrecognised.
 */
inline std::shared_ptr<Thermo> makeThermo(PHARE::initializer::PHAREDict const& dict)
{
    auto const model = getThermoModelFromString(dict["eos"].template to<std::string>());

    switch (model)
    {
        case ThermoModel::ideal_gas:
        {
            double const gamma = dict["heat_capacity_ratio"].template to<double>();
            return std::make_shared<IdealGasThermo>(gamma);
        }
    }

    throw std::runtime_error("makeThermo: unhandled ThermoModel enumerator");
}

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_THERMO_THERMO_FACTORY_HPP
