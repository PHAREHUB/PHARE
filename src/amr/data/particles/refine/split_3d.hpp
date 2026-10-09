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
struct LimePattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<8>>
{
    using Super = SplitPattern<DimConst<3>, RefinedParticlesConst<8>>;

    constexpr LimePattern(float const weight, float const delta)
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
struct PurplePattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<12>>
{
    using Super = SplitPattern<DimConst<3>, RefinedParticlesConst<12>>;

    constexpr PurplePattern(float const weight, float const delta)
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
        Super::deltas_[8]  = {delta, zero, delta};
        Super::deltas_[9]  = {delta, zero, -delta};
        Super::deltas_[10] = {-delta, zero, delta};
        Super::deltas_[11] = {-delta, zero, -delta};
    }
};


template<> // 2 per edge - equidistant from center
struct WhitePattern<DimConst<3>> : SplitPattern<DimConst<3>, RefinedParticlesConst<24>>
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
using SplitPattern_3_1_12_Dispatcher = PatternDispatcher<PurplePattern<DimConst<3>>>;

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
                        LimePattern<DimConst<3>>, PurplePattern<DimConst<3>>>;

template<>
struct Splitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<27>>
    : public ASplitter<DimConst<3>, InterpConst<1>, RefinedParticlesConst<27>>,
      SplitPattern_3_1_27_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_3_1_27_Dispatcher{
              {weight[0]}, {weight[1], delta[0]}, {weight[2], delta[0]}, {weight[3], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta = {1};
    // exact split: tensor product of the 1D {0.5, 0.25} split, ordered as the dispatcher
    // (centre, 6 faces, 8 corners, 12 edges)
    static constexpr std::array<float, 4> weight = {0.125, 0.0625, 0.015625, 0.03125};
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
using SplitPattern_3_2_12_Dispatcher = PatternDispatcher<PurplePattern<DimConst<3>>>;

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
// TODO(@rochSmets) please review and explain: Smets et al. 2021 (CPC, arXiv:2104.10675) Tables 5-6
// give only the tau = 0 and tau = 1 rows for N = 27 at p = 2 and 3. These values are applied to the
// face (Pink), corner (Lime) and edge (Purple) children alike, so the children weights sum
// to 1.538, not 1. Where do these w and delta come from, and what are the tau = 2 and 3 values? Not
// built until then: commented out in res/sim/all.txt (see PHAREHUB/PHARE#1340).
using SplitPattern_3_2_27_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<3>>, PinkPattern<DimConst<3>>,
                        LimePattern<DimConst<3>>, PurplePattern<DimConst<3>>>;

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
using SplitPattern_3_3_12_Dispatcher = PatternDispatcher<PurplePattern<DimConst<3>>>;

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
// TODO(@rochSmets) please review and explain: Smets et al. 2021 (CPC, arXiv:2104.10675) Tables 5-6
// give only the tau = 0 and tau = 1 rows for N = 27 at p = 2 and 3. These values are applied to the
// face (Pink), corner (Lime) and edge (Purple) children alike, so the children weights sum
// to 1.551, not 1. Where do these w and delta come from, and what are the tau = 2 and 3 values? Not
// built until then: commented out in res/sim/all.txt (see PHAREHUB/PHARE#1340).
using SplitPattern_3_3_27_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<3>>, PinkPattern<DimConst<3>>,
                        LimePattern<DimConst<3>>, PurplePattern<DimConst<3>>>;

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
