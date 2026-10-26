#ifndef PHARE_LOHNER_CROSS_DERIVATIVES_TAGGER_STRATEGY_HPP
#define PHARE_LOHNER_CROSS_DERIVATIVES_TAGGER_STRATEGY_HPP

#include "initializer/data_provider.hpp"

#include "per_quantity_tagger_strategy.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>

namespace PHARE::amr
{
/**
 * @brief Löhner's criterion in the form of FLASH (user guide
 * https://flash.rochester.edu/site/flashcode/user_support/flash4_ug_4p62/node59.html, eq. 8.3-8.4),
 * with all cross derivatives, on a stencil of points two cells apart.
 *
 * For a vector quantity, the absolute value of the formula become the L2 norm
 *
 * @tparam Model The physical model
 */
template<typename Model>
class LohnerCrossDerivativesTaggerStrategy
    : public PerQuantityTaggerStrategy<LohnerCrossDerivativesTaggerStrategy<Model>, Model>
{
    using Super = PerQuantityTaggerStrategy<LohnerCrossDerivativesTaggerStrategy<Model>, Model>;

public:
    using recentered_field_type = Super::recentered_field_type;

    static constexpr auto dimension              = Model::dimension; ///< dimension of the model
    static constexpr std::uint32_t stencil_reach = 2; ///< the estimator reads two cells away

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging"
     */
    explicit LohnerCrossDerivativesTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
        , reltol_{cppdict::get_value(dict, "params/reltol", 0.01)}
        , abstol_{cppdict::get_value(dict, "params/abstol", 1e-8)}
    {
    }

    /**
     * @brief The indicator
     *
     * @param components The components of the quantity, viewed at cell centres
     * @param cell The index of the cell in the field arrays
     * @param layout The layout of the patch
     * @return Lohner's indicator with cross derivatives
     */
    double indicator(std::span<recentered_field_type const> components,
                     std::array<std::uint32_t, dimension> const& cell,
                     Model::gridlayout_type const& layout) const
    {
        auto const& meshSize = layout.meshSize();
        auto const shifted   = [&](std::size_t p, int offsetP, std::size_t q, int offsetQ) {
            auto index = cell;
            index[p]   = static_cast<std::uint32_t>(static_cast<int>(index[p]) + offsetP);
            index[q]   = static_cast<std::uint32_t>(static_cast<int>(index[q]) + offsetQ);
            return index;
        };

        double numerator   = 0.;
        double denominator = 0.;
        for (std::size_t p = 0; p < dimension; ++p)
            for (std::size_t q = 0; q < dimension; ++q)
            {
                auto const upperUpIndex   = shifted(p, +1, q, +1);
                auto const lowerUpIndex   = shifted(p, -1, q, +1);
                auto const upperDownIndex = shifted(p, +1, q, -1);
                auto const lowerDownIndex = shifted(p, -1, q, -1);
                double change             = 0.;
                double slopeUp2           = 0.;
                double slopeDown2         = 0.;
                double upperUp2           = 0.;
                double lowerUp2           = 0.;
                double upperDown2         = 0.;
                double lowerDown2         = 0.;
                for (auto const& field : components)
                {
                    auto const upperUp   = field(upperUpIndex);
                    auto const lowerUp   = field(lowerUpIndex);
                    auto const upperDown = field(upperDownIndex);
                    auto const lowerDown = field(lowerDownIndex);
                    auto const slopeUp   = upperUp - lowerUp;
                    auto const slopeDown = upperDown - lowerDown;
                    change += (slopeUp - slopeDown) * (slopeUp - slopeDown);
                    slopeUp2 += slopeUp * slopeUp;
                    slopeDown2 += slopeDown * slopeDown;
                    upperUp2 += upperUp * upperUp;
                    lowerUp2 += lowerUp * lowerUp;
                    upperDown2 += upperDown * upperDown;
                    lowerDown2 += lowerDown * lowerDown;
                }
                auto const filter = std::sqrt(upperUp2) + std::sqrt(lowerUp2)
                                    + std::sqrt(upperDown2) + std::sqrt(lowerDown2);
                auto const weight = 1. / (meshSize[p] * meshSize[q]);
                auto const n      = weight * std::sqrt(change);
                auto const D
                    = weight
                      * (std::sqrt(slopeUp2) + std::sqrt(slopeDown2) + reltol_ * filter + abstol_);
                numerator += n * n;
                denominator += D * D;
            }
        return denominator > 0. ? std::sqrt(numerator / denominator) : 0.;
    }

private:
    double reltol_; ///< weight of the noise filter, relative to the local values
    double abstol_; ///< absolute floor added to the denominator of each pair of directions
};

} // namespace PHARE::amr

#endif
