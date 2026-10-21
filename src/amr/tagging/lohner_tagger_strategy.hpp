#ifndef PHARE_LOHNER_TAGGER_STRATEGY_HPP
#define PHARE_LOHNER_TAGGER_STRATEGY_HPP

#include "initializer/data_provider.hpp"

#include "component_wise_tagger_strategy.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace PHARE::amr
{
template<typename Model>
class LohnerTaggerStrategy : public ComponentWiseTaggerStrategy<LohnerTaggerStrategy<Model>, Model>
{
    using Super = ComponentWiseTaggerStrategy<LohnerTaggerStrategy<Model>, Model>;

public:
    using recentered_field_type                  = Super::recentered_field_type;
    static constexpr auto dimension              = Model::dimension;
    static constexpr std::uint32_t stencil_reach = 1;

    explicit LohnerTaggerStrategy(initializer::PHAREDict const& dict)
        : Super{dict}
        , reltol_{cppdict::get_value(dict, "params/reltol", 0.02)}
        , abstol_{cppdict::get_value(dict, "params/abstol", 1e-30)}
    {
    }

    double indicator(recentered_field_type const& field,
                     std::array<std::uint32_t, dimension> const& cell) const
    {
        auto const a_0 = field(cell);
        double crit    = 0.;
        for (std::size_t d = 0; d < dimension; ++d)
        {
            auto im = cell;
            auto ip = cell;
            im[d] -= 1;
            ip[d] += 1;
            auto const a_m = field(im);
            auto const a_p = field(ip);
            auto const num = std::abs(a_p - 2.0 * a_0 + a_m);
            auto const d1  = std::abs(a_p - a_0) + std::abs(a_0 - a_m);
            auto const flt = std::abs(a_m) + 2.0 * std::abs(a_0) + std::abs(a_p);
            crit           = std::max(crit, num / (d1 + reltol_ * flt + abstol_));
        }
        return crit;
    }

private:
    double reltol_;
    double abstol_;
};

} // namespace PHARE::amr

#endif
