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
    constexpr PatternDispatcher(Patterns&&... patterns)
        : _patterns{patterns...}
        , _nbRefinedPart{(patterns.deltas_.size() + ...)}
    {
    }

    template<size_t dimension, typename Particle, typename RefinedDeltas>
    static void _set_(Particle& particle, RefinedDeltas const& refinedDeltas)
    {
        for (size_t iDim = 0; iDim < dimension; iDim++)
        {
            particle.delta[iDim] += refinedDeltas[iDim];
            float integra = std::floor(particle.delta[iDim]);
            particle.delta[iDim] -= integra;
            particle.iCell[iDim] += static_cast<int32_t>(integra);
        }
    }


    template<size_t dimension, typename Particle, typename Particles>
    void _dispatch(Particle const& particle, Particles& particles, size_t idx) const
    {
        using FineParticle
            = std::conditional_t<std::is_same_v<Particles, std::vector<core::Particle<dimension>>>,
                                 core::Particle<dimension>&,
                                 core::ContiguousParticleView<dimension>>;

        core::apply(_patterns, [&](auto const& pattern) {
            for (size_t rpIndex = 0; rpIndex < pattern.deltas_.size(); rpIndex++)
            {
                FineParticle fineParticle = particles[idx++];
                fineParticle.weight       = particle.weight * pattern.weight_;
                fineParticle.charge       = particle.charge;
                fineParticle.iCell        = particle.iCell;
                fineParticle.delta        = particle.delta;
                fineParticle.v            = particle.v;
                _set_<dimension>(fineParticle, pattern.deltas_[rpIndex]);
            }
        });
    }

    template<size_t dimension, bool OwnedState>
    inline void operator()(core::ContiguousParticleView<dimension> const& coarsePartOnRefinedGrid,
                           core::ContiguousParticles<dimension, OwnedState>& refinedParticles,
                           size_t idx) const
    {
        _dispatch<dimension>(coarsePartOnRefinedGrid, refinedParticles, idx);
    }

    template<size_t dimension>
    inline void operator()(core::Particle<dimension> const& coarsePartOnRefinedGrid,
                           std::vector<core::Particle<dimension>>& refinedParticles) const
    {
        size_t idx = refinedParticles.size();
        for (size_t i = 0; i < _nbRefinedPart; i++)
            refinedParticles.emplace_back();
        _dispatch<dimension>(coarsePartOnRefinedGrid, refinedParticles, idx);
    }

    std::tuple<Patterns...> _patterns{};
    size_t _nbRefinedPart{0};
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
