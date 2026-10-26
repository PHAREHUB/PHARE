#ifndef PHARE_LOHNER_PER_DIRECTION_TAGGER_STRATEGY_HPP
#define PHARE_LOHNER_PER_DIRECTION_TAGGER_STRATEGY_HPP

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
 * @brief Löhner's criterion in the form of PLUTO (Mignone et al. 2012, ApJS 198, 7, eq. 50-51,
 * https://iopscience.iop.org/article/10.1088/0067-0049/198/1/7/pdf).
 *
 * For a vector quantity, the absolute value of the formula become the L2 norm
 *
 * @tparam Model The physical model
 */
template<typename Model>
class LohnerPerDirectionTaggerStrategy
    : public PerQuantityTaggerStrategy<LohnerPerDirectionTaggerStrategy<Model>, Model>
{
    using Super = PerQuantityTaggerStrategy<LohnerPerDirectionTaggerStrategy<Model>, Model>;

public:
    using recentered_field_type = Super::recentered_field_type;

    static constexpr auto dimension              = Model::dimension; ///< dimension of the model
    static constexpr std::uint32_t stencil_reach = 1; ///< the estimator reads the two neighbours

    /**
     * @brief Constructor
     *
     * @param dict config dict at path "simulation/AMR/refinement/tagging"
     */
    explicit LohnerPerDirectionTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
        , reltol_{cppdict::get_value(dict, "params/reltol", 0.01)}
        , abstol_{cppdict::get_value(dict, "params/abstol", 1e-8)}
    {
    }

    /**
     * @brief The indicator of one quantity in one cell
     *
     * @param components The components of the quantity, viewed at cell centres
     * @param cell The index of the cell in the field arrays, the first ghost cell having index 0
     * @param layout unused by this criterion
     * @return indicator value
     */
    double indicator(std::span<recentered_field_type const> components,
                     std::array<std::uint32_t, dimension> const& cell,
                     Model::gridlayout_type const&) const
    {
        double numerator   = 0.;
        double denominator = 0.;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto im = cell;
            auto ip = cell;
            im[d] -= 1;
            ip[d] += 1;
            double change = 0.;
            double slopeM = 0.;
            double slopeP = 0.;
            double normM  = 0.;
            double norm0  = 0.;
            double normP  = 0.;
            for (auto const& field : components)
            {
                auto const a_m = field(im);
                auto const a_0 = field(cell);
                auto const a_p = field(ip);
                change += (a_p - 2.0 * a_0 + a_m) * (a_p - 2.0 * a_0 + a_m);
                slopeM += (a_0 - a_m) * (a_0 - a_m);
                slopeP += (a_p - a_0) * (a_p - a_0);
                normM += a_m * a_m;
                norm0 += a_0 * a_0;
                normP += a_p * a_p;
            }
            auto const D
                = std::sqrt(slopeP) + std::sqrt(slopeM)
                  + reltol_ * (std::sqrt(normM) + 2.0 * std::sqrt(norm0) + std::sqrt(normP))
                  + abstol_;
            numerator += change;
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
