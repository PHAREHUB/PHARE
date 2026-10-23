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
 * @brief TOFILL
 *
 * @tparam Model TOFILL
 */
template<typename Model>
class DefaultTaggerStrategy
    : public ComponentWiseTaggerStrategy<DefaultTaggerStrategy<Model>, Model>
{
    using Super = ComponentWiseTaggerStrategy<DefaultTaggerStrategy<Model>, Model>;

public:
    /// TOFILL
    using recentered_field_type                  = Super::recentered_field_type;
    static constexpr auto dimension              = Model::dimension;       ///< TOFILL
    static constexpr std::uint32_t stencil_reach = dimension == 1 ? 3 : 2; ///< TOFILL

    /**
     * @brief TOFILL
     *
     * @param dict TOFILL
     */
    explicit DefaultTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
    {
    }

    /**
     * @brief TOFILL
     *
     * @param F TOFILL
     * @param cell TOFILL
     * @param layout unused by this criterion
     * @return TOFILL
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
     * @brief TOFILL
     *
     * @param indicators TOFILL
     * @return TOFILL
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
