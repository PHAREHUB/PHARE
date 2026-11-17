#ifndef PHARE_LOHNER_PARAMETERS_HPP
#define PHARE_LOHNER_PARAMETERS_HPP

#include "initializer/data_provider.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace PHARE::amr
{
/**
 * @brief The per-quantity parameters of both forms of Löhner's criterion.
 */
struct LohnerParameters
{
    double reltol = 0.01; ///< weight of the noise filter, relative to the local values
    double refMagnitude
        = 1e-8; ///< reference magnitude added to the denominator, in the units of the quantity

    /**
     * @brief Read the parameters of one tagged quantity.
     *
     * @param quantity config dict of the quantity, at path "simulation/AMR/refinement/tagging/Q<i>"
     * @param name The name of the quantity, for error messages
     * @return the parameters, with their defaults where absent
     * @throws std::runtime_error if reltol or ref_magnitude is negative or not finite
     */
    static LohnerParameters read(initializer::PHAREDict const& quantity, std::string const& name)
    {
        LohnerParameters constexpr defaults{};
        LohnerParameters const parameters{
            cppdict::get_value(quantity, "reltol", defaults.reltol),
            cppdict::get_value(quantity, "ref_magnitude", defaults.refMagnitude)};
        auto const check = [&](std::string const& key, double value) {
            if (!(std::isfinite(value) and value >= 0.))
                throw std::runtime_error("tagging quantity '" + name + "': " + key
                                         + " must be finite and >= 0, got "
                                         + std::to_string(value));
        };
        check("reltol", parameters.reltol);
        check("ref_magnitude", parameters.refMagnitude);
        return parameters;
    }
};

} // namespace PHARE::amr

#endif
