#ifndef PHARE_LOHNER_PER_DIRECTION_TAGGER_STRATEGY_HPP
#define PHARE_LOHNER_PER_DIRECTION_TAGGER_STRATEGY_HPP

#include "initializer/data_provider.hpp"

#include "component_wise_tagger_strategy.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace PHARE::amr
{
/**
 * @brief Löhner's error estimator in the form of PLUTO (Mignone et al. 2012, ApJS 198, 7,
 * eq. 50-51), without cross derivatives.
 *
 * Along each direction d, the second difference n_d = |a+ - 2 a0 + a-| is compared with
 * D_d = |a+ - a0| + |a0 - a-| + reltol (|a+| + 2 |a0| + |a-|) + abstol, and the directions are
 * combined in a single ratio, sqrt(sum_d n_d^2 / sum_d D_d^2), which lies in [0, 1/(1 + reltol)).
 *
 * @tparam Model The physical model
 */
template<typename Model>
class LohnerPerDirectionTaggerStrategy
    : public ComponentWiseTaggerStrategy<LohnerPerDirectionTaggerStrategy<Model>, Model>
{
    using Super = ComponentWiseTaggerStrategy<LohnerPerDirectionTaggerStrategy<Model>, Model>;

public:
    /// cell-centred view of one component of a tagged quantity
    using recentered_field_type                  = Super::recentered_field_type;
    static constexpr auto dimension              = Model::dimension; ///< dimension of the model
    static constexpr std::uint32_t stencil_reach = 1; ///< the estimator reads the two neighbours

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging", from which
     * "params/reltol" (default 0.01) and "params/abstol" (default 1e-8) are read
     */
    explicit LohnerPerDirectionTaggerStrategy(initializer::PHAREDict const& dict)
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
     * @param layout unused by this criterion
     * @return the ratio of summed squares, 0 when the denominator vanishes
     */
    double indicator(recentered_field_type const& field,
                     std::array<std::uint32_t, dimension> const& cell,
                     Model::gridlayout_type const&) const
    {
        auto const a_0     = field(cell);
        double numerator   = 0.;
        double denominator = 0.;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto im = cell;
            auto ip = cell;
            im[d] -= 1;
            ip[d] += 1;
            auto const a_m = field(im);
            auto const a_p = field(ip);
            auto const n   = std::abs(a_p - 2.0 * a_0 + a_m);
            auto const D   = std::abs(a_p - a_0) + std::abs(a_0 - a_m)
                           + reltol_ * (std::abs(a_m) + 2.0 * std::abs(a_0) + std::abs(a_p))
                           + abstol_;
            numerator += n * n;
            denominator += D * D;
        }
        return denominator > 0. ? std::sqrt(numerator / denominator) : 0.;
    }

private:
    double reltol_; ///< weight of the noise filter, relative to the local values
    double abstol_; ///< absolute floor added to the denominator of each direction
};

} // namespace PHARE::amr

#endif
