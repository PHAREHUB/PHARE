#ifndef DEFAULT_TAGGER_STRATEGY_H
#define DEFAULT_TAGGER_STRATEGY_H

#include "initializer/data_provider.hpp"

#include "component_wise_tagger_strategy.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <algorithm>

namespace PHARE::amr
{
/**
 * @brief PHARE's default refinement criterion.
 *
 * @tparam Model The physical model
 */
template<typename Model>
class DefaultTaggerStrategy
    : public ComponentWiseTaggerStrategy<DefaultTaggerStrategy<Model>, Model>
{
    using Super = ComponentWiseTaggerStrategy<DefaultTaggerStrategy<Model>, Model>;

public:
    /// cell-centred view of one component of a tagged quantity
    using recentered_field_type = Super::recentered_field_type;

    static constexpr auto dimension = Model::dimension; ///< dimension of the model
    /// the indicator reads F(i-2) to F(i+3) in 1D, F(i) to F(i+2) in 2D and 3D
    static constexpr std::uint32_t stencil_reach = dimension == 1 ? 3 : 2;

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging"
     */
    explicit DefaultTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    /**
     * @brief The indicator of one component in one cell
     *
     * @param F The component, viewed at cell centres
     * @param cell The index of the cell in the field arrays, the first ghost cell having index 0
     * @param layout unused by this criterion
     * @return the indicator value at the cell
     */
    double indicator(recentered_field_type const& F,
                     std::array<std::uint32_t, dimension> const& cell,
                     Model::gridlayout_type const&) const
    {
        if constexpr (dimension == 1)
        {
            auto const ix     = cell[0];
            auto const Favg   = 0.2 * (F(ix - 2) + F(ix - 1) + F(ix) + F(ix + 1) + F(ix + 2));
            auto const Favgp1 = 0.2 * (F(ix - 1) + F(ix) + F(ix + 1) + F(ix + 2) + F(ix + 3));
            return std::abs(Favgp1 - Favg) / (1 + std::abs(Favg));
        }
        else
        {
            auto const f0 = F(cell);
            double crit   = 0.;
            for (std::size_t d = 0; d < dimension; ++d)
            {
                auto next  = cell;
                auto next2 = cell;
                next[d] += 1;
                next2[d] += 2;
                auto const delta_2 = std::abs(F(next2) - f0);
                auto const delta   = std::abs(F(next) - f0);
                crit               = std::max(crit, delta_2 / (1 + delta));
            }
            return crit;
        }
    }

    /**
     * @brief Combines the indicators of the components of a quantity, in 1D only
     *
     * @param indicators The indicator of each component of the quantity
     * @return their L2 norm, compared with the quantity's threshold
     */
    static double combine(std::span<double const> indicators)
        requires(dimension == 1)
    {
        double sum = 0.;
        for (auto const criter : indicators)
            sum += criter * criter;
        return std::sqrt(sum);
    }
};

} // namespace PHARE::amr

#endif
