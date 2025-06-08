/*
Splitting reference material can be found @
  https://github.com/PHAREHUB/PHARE/wiki/SplitPattern
*/

#ifndef PHARE_SPLIT_3D_HPP
#define PHARE_SPLIT_3D_HPP

#include <array>
#include <cstddef>
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"
#include "splitter.hpp"

namespace PHARE::amr
{
using namespace PHARE::core;

/**************************************************************************/
template<> // 1 per face - centered
struct PinkPattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<6>>
{
    using Super = SplitPattern<DimConst<3>, RefinedParticlesConst<6>>;

    constexpr PinkPattern(float const weight, float const delta)
        : Super{weight}
    {
        constexpr float zero = 0;

        for (std::size_t i = 0; i < 2; i++)
        {
            std::size_t offset = i * 3;
            float sign         = i % 2 ? -1 : 1;
            auto mode          = delta * sign;

            Super::deltas_[0 + offset] = {mode, zero, zero};
            Super::deltas_[1 + offset] = {zero, zero, mode};
            Super::deltas_[2 + offset] = {zero, mode, zero};
        }
    }
};

template<> // 1 per corner
struct PurplePattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<8>>
{
    using Super = SplitPattern<DimConst<3>, RefinedParticlesConst<8>>;

    constexpr PurplePattern(float const weight, float const delta)
        : Super{weight}
    {
        for (std::size_t i = 0; i < 2; i++)
        {
            std::size_t offset = i * 4;
            float sign         = i % 2 ? -1 : 1;
            auto mode          = delta * sign;

            Super::deltas_[0 + offset] = {mode, mode, mode};
            Super::deltas_[1 + offset] = {mode, mode, -mode};
            Super::deltas_[2 + offset] = {mode, -mode, mode};
            Super::deltas_[3 + offset] = {mode, -mode, -mode};
        }
    }
};


template<> // 1 per edge - centered
struct LimePattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<12>>
{
    using Super = SplitPattern<DimConst<3>, RefinedParticlesConst<12>>;

    constexpr LimePattern(float const weight, float const delta)
        : Super{weight}
    {
        constexpr float zero = 0;

        auto addSquare = [&delta, this](size_t offset, float y) {
            Super::deltas_[0 + offset] = {zero, y, delta};
            Super::deltas_[1 + offset] = {zero, y, -delta};
            Super::deltas_[2 + offset] = {delta, y, zero};
            Super::deltas_[3 + offset] = {-delta, y, zero};
        };

        addSquare(0, delta);  // top
        addSquare(4, -delta); // bottom
        addSquare(8, 0);      // middle
    }
};


template<> // 2 per edge - equidistant from center
class WhitePattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<24>>
{
    using Super = SplitPattern<DimConst<3>, RefinedParticlesConst<24>>;

    template<typename Delta>
    constexpr WhitePattern(float const weight, Delta const& delta)
        : Super{weight}
    {
        auto addSquare = [&delta, this](size_t offset, float x, float y, float z) {
            Super::deltas_[0 + offset] = {x, y, z};
            Super::deltas_[1 + offset] = {x, -y, z};
            Super::deltas_[2 + offset] = {-x, y, z};
            Super::deltas_[3 + offset] = {-x, -y, z};
        };

        auto addSquares = [addSquare, &delta, this](size_t offset, float x, float y, float z) {
            addSquare(offset, x, y, z);
            addSquare(offset + 4, x, y, -z);
        };

        addSquares(0, delta[0], delta[0], delta[1]);
        addSquares(8, delta[0], delta[1], delta[0]);
        addSquares(16, delta[1], delta[0], delta[0]);
    }
};


/**************************************************************************/
/****************************** INTERP == 1 *******************************/
/**************************************************************************/
using SplitPattern_3_1_6_Dispatcher = PatternDispatcher<PinkPattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<6>>
    : public ASplitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<6>>,
      SplitPattern_3_1_6_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_1_6_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {.966431};
    static constexpr std::array<float, 1> weight = {.166666};
};


/**************************************************************************/
using SplitPattern_3_1_12_Dispatcher = PatternDispatcher<LimePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<12>>
    : public ASplitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<12>>,
      SplitPattern_3_1_12_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_1_12_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {.74823};
    static constexpr std::array<float, 1> weight = {.083333};
};


/**************************************************************************/
using SplitPattern_3_1_27_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<3>>, PinkPattern<DimConst<3>>,
                        PurplePattern<DimConst<3>>, LimePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<27>>
    : public ASplitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<27>>,
      SplitPattern_3_1_27_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_1_27_Dispatcher{
            {weight[0]}, {weight[1], delta[0]}, {weight[1], delta[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1};
    static constexpr std::array<float, 2> weight = {0.125, 0.0625};
};



/**************************************************************************/
/****************************** INTERP == 2 *******************************/
/**************************************************************************/
using SplitPattern_3_2_6_Dispatcher = PatternDispatcher<PinkPattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<2>, RefinedParticlesConst<6>>
    : public ASplitter<DimConst<3>, InterpConst<2>, RefinedParticlesConst<6>>,
      SplitPattern_3_2_6_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_2_6_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.149658};
    static constexpr std::array<float, 1> weight = {.166666};
};


/**************************************************************************/
using SplitPattern_3_2_12_Dispatcher = PatternDispatcher<LimePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<2>, RefinedParticlesConst<12>>
    : public ASplitter<DimConst<3>, InterpConst<2>, RefinedParticlesConst<12>>,
      SplitPattern_3_2_12_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_2_12_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {.888184};
    static constexpr std::array<float, 1> weight = {.083333};
};


/**************************************************************************/
using SplitPattern_3_2_27_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<3>>, PinkPattern<DimConst<3>>,
                        PurplePattern<DimConst<3>>, LimePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<2>, RefinedParticlesConst<27>>
    : public ASplitter<DimConst<3>, InterpConst<2>, RefinedParticlesConst<27>>,
      SplitPattern_3_2_27_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_2_27_Dispatcher{
            {weight[0]}, {weight[1], delta[0]}, {weight[1], delta[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.111333};
    static constexpr std::array<float, 2> weight = {.099995, .055301};
};


/**************************************************************************/
/****************************** INTERP == 3 *******************************/
/**************************************************************************/
using SplitPattern_3_3_6_Dispatcher = PatternDispatcher<PinkPattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<3>, RefinedParticlesConst<6>>
    : public ASplitter<DimConst<3>, InterpConst<3>, RefinedParticlesConst<6>>,
      SplitPattern_3_3_6_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_3_6_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.312622};
    static constexpr std::array<float, 1> weight = {.166666};
};


/**************************************************************************/
using SplitPattern_3_3_12_Dispatcher = PatternDispatcher<LimePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<3>, RefinedParticlesConst<12>>
    : public ASplitter<DimConst<3>, InterpConst<3>, RefinedParticlesConst<12>>,
      SplitPattern_3_3_12_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_3_12_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.012756};
    static constexpr std::array<float, 1> weight = {.083333};
};


/**************************************************************************/
using SplitPattern_3_3_27_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<3>>, PinkPattern<DimConst<3>>,
                        PurplePattern<DimConst<3>>, LimePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<3>, RefinedParticlesConst<27>>
    : public ASplitter<DimConst<3>, InterpConst<3>, RefinedParticlesConst<27>>,
      SplitPattern_3_3_27_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_3_27_Dispatcher{
            {weight[0]}, {weight[1], delta[0]}, {weight[1], delta[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.276815};
    static constexpr std::array<float, 2> weight = {.104047, .05564};
};


/**************************************************************************/


} // namespace PHARE::amr


#endif /* PHARE_SPLIT_3D_HPP */
