#ifndef PHARE_SPLIT_H
#define PHARE_SPLIT_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "core/utilities/types.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/particles/particle.h"
#include "core/utilities/box/box.h"


namespace PHARE::amr
{
// see meta_utilities.h for list of declared permutations.

template<typename _dimension, typename _nbRefinedPart>
struct SplitPattern
{
    constexpr SplitPattern(float weight)
        : weight_{weight}
    {
    }

    float weight_;
    std::array<core::Point<float, _dimension{}()>, _nbRefinedPart{}()> deltas_{};
};

template<typename... Patterns>
struct PatternDispatcher
{
    constexpr PatternDispatcher(Patterns&&... patterns_)
        : patterns{patterns_...}
    {
    }

    template<std::size_t dimension>
    inline void operator()(core::Particle<dimension> const& coarsePartOnRefinedGrid,
                           std::vector<core::Particle<dimension>>& refinedParticles) const
    {
        auto get = [](size_t iDim, auto& icell, auto& delta, auto& refinedDeltas) {
            delta[iDim] += refinedDeltas[iDim];
            float integra = std::floor(delta[iDim]);
            delta[iDim] -= integra;
            icell[iDim] += static_cast<int32_t>(integra);
        };

        core::apply(patterns, [&](auto const& pattern) {
            for (size_t rpIndex = 0; rpIndex < pattern.deltas_.size(); rpIndex++)
            {
                auto icell = coarsePartOnRefinedGrid.iCell;
                auto delta = coarsePartOnRefinedGrid.delta;

                for (size_t i = 0; i < dimension; i++)
                    get(i, icell, delta, pattern.deltas_[rpIndex]);

                float weight = coarsePartOnRefinedGrid.weight * pattern.weight_;
                refinedParticles.push_back({weight, coarsePartOnRefinedGrid.charge, icell, delta,
                                            coarsePartOnRefinedGrid.v});
            }
        });
    }

    std::tuple<Patterns...> patterns{};
};


template<typename _dimension, typename _interp_order, typename _nbRefinedPart>
struct ASplitter
{
    static constexpr auto dimension     = _dimension{}();
    static constexpr auto interp_order  = _interp_order{}();
    static constexpr auto nbRefinedPart = _nbRefinedPart{}();

    static constexpr int maxCellDistanceFromSplit()
    {
        constexpr auto particleSize = interp_order + 1;
        return std::ceil(particleSize * 0.5);
    }

    constexpr ASplitter() {}
};

template<typename dimension, typename interp_order, typename nbRefinedPart>
class Splitter : public ASplitter<dimension, interp_order, nbRefinedPart>
{
    Splitter() = delete; // Unspecialized template class, never to be instantiated
};

template<typename dim>
struct BlackDispatcher : SplitPattern<dim, core::RefinedParticlesConst<1>>
{
    using Super = SplitPattern<dim, core::RefinedParticlesConst<1>>;
    constexpr BlackDispatcher(float const weight)
        : Super{weight}
    {
    }
};

template<typename dim>
struct PurpleDispatcher
{
};
template<typename dim>
struct BrownDispatcher
{
};
template<typename dim>
struct PinkDispatcher
{
};

} // namespace PHARE::amr

#include "split_1d.h"
#include "split_2d.h"


#endif // endif SPLIT_H
