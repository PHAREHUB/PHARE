#ifndef PHARE_LOHNER_CROSS_DERIVATIVES_TAGGER_STRATEGY_HPP
#define PHARE_LOHNER_CROSS_DERIVATIVES_TAGGER_STRATEGY_HPP

#include "initializer/data_provider.hpp"

#include "component_wise_tagger_strategy.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace PHARE::amr
{
/**
 * @brief Löhner's error estimator in the form of FLASH (user guide §8.6.3.1, eq. 8.3-8.4), with
 * all cross derivatives, on a stencil of points two cells apart.
 *
 * With the slope over two cells S_p(y) = u(y + e_p) - u(y - e_p), every ordered pair of directions
 * (p, q) contributes n_pq = |S_p(x + e_q) - S_p(x - e_q)| and
 * D_pq = |S_p(x + e_q)| + |S_p(x - e_q)| + reltol (sum of |u| over the four points used) + abstol,
 * both weighted by 1/(dx_p dx_q). The indicator sqrt(sum (w n)^2 / sum (w D)^2) lies in
 * [0, 1/(1 + reltol)), and only depends on the ratios between the mesh sizes.
 *
 * @tparam Model The physical model
 */
template<typename Model>
class LohnerCrossDerivativesTaggerStrategy
    : public ComponentWiseTaggerStrategy<LohnerCrossDerivativesTaggerStrategy<Model>, Model>
{
    using Super = ComponentWiseTaggerStrategy<LohnerCrossDerivativesTaggerStrategy<Model>, Model>;

public:
    /// cell-centred view of one component of a tagged quantity
    using recentered_field_type                  = Super::recentered_field_type;
    static constexpr auto dimension              = Model::dimension; ///< dimension of the model
    static constexpr std::uint32_t stencil_reach = 2; ///< the estimator reads two cells away

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging", from which
     * "params/reltol" (default 0.01) and "params/abstol" (default 1e-8) are read
     */
    explicit LohnerCrossDerivativesTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
        , reltol_{cppdict::get_value(dict, "params/reltol", 0.01)}
        , abstol_{cppdict::get_value(dict, "params/abstol", 1e-8)}
    {
    }

    /**
     * @brief The estimator of one component in one cell
     *
     * @param field The component, viewed at cell centres
     * @param cell The index of the cell in the field arrays
     * @param layout The layout of the patch, whose mesh sizes weight each pair of directions
     * @return the ratio of weighted summed squares, 0 when the denominator vanishes
     */
    double indicator(recentered_field_type const& field,
                     std::array<std::uint32_t, dimension> const& cell,
                     Model::gridlayout_type const& layout) const
    {
        auto const& meshSize = layout.meshSize();
        auto const value     = [&](std::size_t p, int offsetP, std::size_t q, int offsetQ) {
            auto index = cell;
            index[p]   = static_cast<std::uint32_t>(static_cast<int>(index[p]) + offsetP);
            index[q]   = static_cast<std::uint32_t>(static_cast<int>(index[q]) + offsetQ);
            return field(index);
        };

        double numerator   = 0.;
        double denominator = 0.;
        for (std::size_t p = 0; p < dimension; ++p)
            for (std::size_t q = 0; q < dimension; ++q)
            {
                auto const upperUp   = value(p, +1, q, +1);
                auto const lowerUp   = value(p, -1, q, +1);
                auto const upperDown = value(p, +1, q, -1);
                auto const lowerDown = value(p, -1, q, -1);
                auto const slopeUp   = upperUp - lowerUp;
                auto const slopeDown = upperDown - lowerDown;
                auto const filter    = std::abs(upperUp) + std::abs(lowerUp) + std::abs(upperDown)
                                    + std::abs(lowerDown);
                auto const weight = 1. / (meshSize[p] * meshSize[q]);
                auto const n      = weight * std::abs(slopeUp - slopeDown);
                auto const D
                    = weight
                      * (std::abs(slopeUp) + std::abs(slopeDown) + reltol_ * filter + abstol_);
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
