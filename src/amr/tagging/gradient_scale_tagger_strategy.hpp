#ifndef PHARE_GRADIENT_SCALE_TAGGER_STRATEGY_HPP
#define PHARE_GRADIENT_SCALE_TAGGER_STRATEGY_HPP

#include "initializer/data_provider.hpp"

#include "amr/tagging/per_quantity_tagger_strategy.hpp"

#include "dict.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>

namespace PHARE::amr
{
/**
 * @brief The per-quantity parameters of the gradient-scale criterion.
 */
struct GradientScaleParameters
{
    double refMagnitude; ///< added to the magnitude of the quantity, in its units

    /**
     * @brief Read the parameters of one tagged quantity.
     *
     * @param quantity config dict of the quantity, at path "simulation/AMR/refinement/tagging/Q<i>"
     * @param name The name of the quantity, for error messages
     * @return the parameters
     * @throws std::runtime_error if ref_magnitude is absent, not positive or not finite
     */
    static GradientScaleParameters read(initializer::PHAREDict const& quantity,
                                        std::string const& name)
    {
        if (!quantity.contains("ref_magnitude"))
            throw std::runtime_error("tagging quantity '" + name
                                     + "': gradient-scale needs a ref_magnitude");
        GradientScaleParameters const parameters{quantity["ref_magnitude"].template to<double>()};
        if (!(std::isfinite(parameters.refMagnitude) and parameters.refMagnitude > 0.))
            throw std::runtime_error("tagging quantity '" + name
                                     + "': ref_magnitude must be finite and > 0, got "
                                     + std::to_string(parameters.refMagnitude));
        return parameters;
    }
};

/**
 * @brief Cell size over the gradient length of a quantity, dx |grad V|_F / (|V| + ref_magnitude).
 *
 * With V_j the components of the quantity (a single one for a scalar), the indicator is
 * sqrt(sum_d sum_j ((V_j(c + e_d) - V_j(c - e_d)) / 2)^2) / (sqrt(sum_j V_j(c)^2) + ref_magnitude),
 * the undivided centred differences carrying the cell size of each direction. On square cells, this
 * is dx |grad V|_F / (|V| + ref_magnitude), with |grad V|_F^2 = sum_ij (d_i V_j)^2.
 *
 * @tparam Model The physical model
 */
template<typename Model>
class GradientScaleTaggerStrategy
    : public PerQuantityTaggerStrategy<GradientScaleTaggerStrategy<Model>, Model,
                                       GradientScaleParameters>
{
    using Super = PerQuantityTaggerStrategy<GradientScaleTaggerStrategy<Model>, Model,
                                            GradientScaleParameters>;

public:
    using recentered_field_type = Super::recentered_field_type;

    static constexpr auto dimension              = Model::dimension; ///< dimension of the model
    static constexpr std::uint32_t stencil_reach = 1; ///< the estimator reads the two neighbours

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging"
     */
    explicit GradientScaleTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    /**
     * @brief The indicator of one quantity in one cell
     *
     * @param components The components of the quantity, viewed at cell centres
     * @param cell The index of the cell in the field arrays, the first ghost cell having index 0
     * @param layout unused by this criterion
     * @param parameters The quantity's ref_magnitude
     * @return indicator value
     */
    double indicator(std::span<recentered_field_type const> components,
                     std::array<std::uint32_t, dimension> const& cell,
                     Model::gridlayout_type const&, GradientScaleParameters const& parameters) const
    {
        double gradient2 = 0.;
        double norm2     = 0.;
        for (auto const& field : components)
        {
            auto const value = field(cell);
            norm2 += value * value;
            for (std::size_t d = 0; d < dimension; ++d)
            {
                auto im = cell;
                auto ip = cell;
                im[d] -= 1;
                ip[d] += 1;
                auto const halfDifference = 0.5 * (field(ip) - field(im));
                gradient2 += halfDifference * halfDifference;
            }
        }
        return std::sqrt(gradient2) / (std::sqrt(norm2) + parameters.refMagnitude);
    }
};

} // namespace PHARE::amr

#endif
