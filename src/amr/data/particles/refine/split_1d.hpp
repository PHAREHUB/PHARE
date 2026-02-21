/*
Splitting reference material can be found @
  https://github.com/PHAREHUB/PHARE/wiki/SplitPattern
*/

#ifndef PHARE_SPLIT_1D_HPP
#define PHARE_SPLIT_1D_HPP

#include <array>
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"
#include "splitter.hpp"

namespace PHARE::amr
{
using namespace PHARE::core;

/**************************************************************************/
template<>
struct PinkPattern<DimConst<1>> : SplitPattern<DimConst<1>, RefinedParticlesConst<2>>
{
    using Super = SplitPattern<DimConst<1>, RefinedParticlesConst<2>>;

    constexpr PinkPattern(float const weight, float const delta)
        : Super{weight}
    {
        Super::deltas_[0] = {delta};
        Super::deltas_[1] = {-delta};
    }
};


/**************************************************************************/
using SplitPattern_1_1_2_Dispatcher = PatternDispatcher<PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<1>, RefinedParticlesConst<2>>
    : public ASplitter<DimConst<1>, InterpConst<1>, RefinedParticlesConst<2>>,
      SplitPattern_1_1_2_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_1_2_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.551569};
    static constexpr std::array<float, 1> weight = {0.5};
};


/**************************************************************************/
using SplitPattern_1_1_3_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<1>>, PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<1>, RefinedParticlesConst<3>>
    : public ASplitter<DimConst<1>, InterpConst<1>, RefinedParticlesConst<3>>,
      SplitPattern_1_1_3_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_1_3_Dispatcher{{weight[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1};
    static constexpr std::array<float, 2> weight = {0.5, 0.25};
};


/**************************************************************************/
using SplitPattern_1_2_2_Dispatcher = PatternDispatcher<PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<2>, RefinedParticlesConst<2>>
    : public ASplitter<DimConst<1>, InterpConst<2>, RefinedParticlesConst<2>>,
      SplitPattern_1_2_2_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_2_2_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.663959};
    static constexpr std::array<float, 1> weight = {0.5};
};


/**************************************************************************/
using SplitPattern_1_2_3_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<1>>, PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<2>, RefinedParticlesConst<3>>
    : public ASplitter<DimConst<1>, InterpConst<2>, RefinedParticlesConst<3>>,
      SplitPattern_1_2_3_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_2_3_Dispatcher{{weight[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.112033};
    static constexpr std::array<float, 2> weight = {0.468137, 0.265931};
};


/**************************************************************************/
using SplitPattern_1_2_4_Dispatcher
    = PatternDispatcher<PinkPattern<DimConst<1>>, PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<2>, RefinedParticlesConst<4>>
    : public ASplitter<DimConst<1>, InterpConst<2>, RefinedParticlesConst<4>>,
      SplitPattern_1_2_4_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_2_4_Dispatcher{{weight[0], delta[0]}, {weight[1], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {0.5, 1.5};
    static constexpr std::array<float, 2> weight = {0.375, 0.125};
};


/**************************************************************************/
using SplitPattern_1_3_2_Dispatcher = PatternDispatcher<PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<2>>
    : public ASplitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<2>>,
      SplitPattern_1_3_2_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_2_2_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.752399};
    static constexpr std::array<float, 1> weight = {0.5};
};


/**************************************************************************/
using SplitPattern_1_3_3_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<1>>, PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<3>>
    : public ASplitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<3>>,
      SplitPattern_1_3_3_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_3_3_Dispatcher{{weight[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.275922};
    static constexpr std::array<float, 2> weight = {0.473943, 0.263028};
};


/**************************************************************************/
using SplitPattern_1_3_4_Dispatcher
    = PatternDispatcher<PinkPattern<DimConst<1>>, PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<4>>
    : public ASplitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<4>>,
      SplitPattern_1_3_4_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_3_4_Dispatcher{{weight[0], delta[0]}, {weight[1], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {0.542949, 1.664886};
    static constexpr std::array<float, 2> weight = {0.364766, 0.135234};
};


/**************************************************************************/
using SplitPattern_1_3_5_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<1>>, PinkPattern<DimConst<1>>,
                        PinkPattern<DimConst<1>>>;

template<>
struct Splitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<5>>
    : public ASplitter<DimConst<1>, InterpConst<3>, RefinedParticlesConst<5>>,
      SplitPattern_1_3_5_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_1_3_5_Dispatcher{{weight[0]}, {weight[1], delta[0]}, {weight[2], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {1, 2};
    static constexpr std::array<float, 3> weight = {0.375, 0.25, 0.0625};
};


} // namespace PHARE::amr
#endif /*PHARE_SPLIT_1D_H*/
