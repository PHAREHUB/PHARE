/*
Splitting reference material can be found @
  https://github.com/PHAREHUB/PHARE/wiki/SplitPattern
*/

#ifndef PHARE_SPLIT_2D_HPP
#define PHARE_SPLIT_2D_HPP

#include <array>
#include <cstddef>
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"
#include "splitter.hpp"

namespace PHARE::amr
{
using namespace PHARE::core;

/**************************************************************************/
template<>
struct PurplePattern<DimConst<2>> : SplitPattern<DimConst<2>, RefinedParticlesConst<4>>
{
    using Super = SplitPattern<DimConst<2>, RefinedParticlesConst<4>>;

    constexpr PurplePattern(float const weight, float const delta)
        : Super{weight}
    {
        Super::deltas_[0] = {-delta, -delta};
        Super::deltas_[1] = {-delta, +delta};
        Super::deltas_[2] = {+delta, -delta};
        Super::deltas_[3] = {+delta, +delta};
    }
};


template<>
struct BrownPattern<DimConst<2>> : SplitPattern<DimConst<2>, RefinedParticlesConst<8>>
{
    using Super = SplitPattern<DimConst<2>, RefinedParticlesConst<8>>;

    template<typename Delta>
    constexpr BrownPattern(float const weight, Delta const& delta)
        : Super{weight}
    {
        for (std::size_t i = 0; i < 2; i++)
        {
            std::size_t offset         = (i * 4);
            float sign                 = i % 2 ? -1 : 1;
            Super::deltas_[0 + offset] = {-delta[0] * sign, delta[1] * sign};
            Super::deltas_[1 + offset] = {-delta[1] * sign, delta[0] * sign};
            Super::deltas_[2 + offset] = {delta[0] * sign, delta[1] * sign};
            Super::deltas_[3 + offset] = {delta[1] * sign, delta[0] * sign};
        }
    }
};

template<>
struct PinkPattern<DimConst<2>> : SplitPattern<DimConst<2>, RefinedParticlesConst<4>>
{
    using Super = SplitPattern<DimConst<2>, RefinedParticlesConst<4>>;

    constexpr PinkPattern(float const weight, float const delta)
        : Super{weight}
    {
        Super::deltas_[0] = {0.0f, -delta};
        Super::deltas_[1] = {-delta, 0.0f};
        Super::deltas_[2] = {+delta, 0.0f};
        Super::deltas_[3] = {0.0f, delta};
    }
};


/**************************************************************************/
using SplitPattern_2_1_4_Dispatcher = PatternDispatcher<PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<4>>
    : public ASplitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<4>>,
      SplitPattern_2_1_4_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_1_4_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.571783};
    static constexpr std::array<float, 1> weight = {0.25};
};


/**************************************************************************/
using SplitPattern_2_1_5_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<5>>
    : public ASplitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<5>>,
      SplitPattern_2_1_5_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_1_5_Dispatcher{{weight[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.721835};
    static constexpr std::array<float, 2> weight = {0.239863, 0.190034};
};


/**************************************************************************/
using SplitPattern_2_1_8_Dispatcher
    = PatternDispatcher<PurplePattern<DimConst<2>>, PinkPattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<8>>
    : public ASplitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<8>>,
      SplitPattern_2_1_8_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_1_8_Dispatcher{{1, 2}, {2, 1}}
    {
    }

    static constexpr std::array<float, 2> delta  = {0.700909, 1.05786};
    static constexpr std::array<float, 2> weight = {0.179488, 0.070512};
};


/**************************************************************************/
using SplitPattern_2_1_9_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PinkPattern<DimConst<2>>,
                        PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<9>>
    : public ASplitter<DimConst<2>, InterpConst<1>, RefinedParticlesConst<9>>,
      SplitPattern_2_1_9_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_1_9_Dispatcher{{weight[0]}, {weight[1], delta[0]}, {weight[2], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1};
    static constexpr std::array<float, 3> weight = {0.25, 0.125, 0.01625};
};


/**************************************************************************/
using SplitPattern_2_2_4_Dispatcher = PatternDispatcher<PinkPattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<4>>
    : public ASplitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<4>>,
      SplitPattern_2_2_4_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_2_4_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.683734};
    static constexpr std::array<float, 1> weight = {0.25};
};


/**************************************************************************/
using SplitPattern_2_2_5_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PinkPattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<5>>
    : public ASplitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<5>>,
      SplitPattern_2_2_5_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_2_5_Dispatcher{{weight[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.203227};
    static constexpr std::array<float, 2> weight = {0.239166, 0.190209};
};


/**************************************************************************/
using SplitPattern_2_2_8_Dispatcher
    = PatternDispatcher<PinkPattern<DimConst<2>>, PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<8>>
    : public ASplitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<8>>,
      SplitPattern_2_2_8_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_2_8_Dispatcher{{weight[0], delta[0]}, {weight[1], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {0.828428, 1.236701};
    static constexpr std::array<float, 2> weight = {0.178624, 0.071376};
};


/**************************************************************************/
using SplitPattern_2_2_9_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PinkPattern<DimConst<2>>,
                        PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<9>>
    : public ASplitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<9>>,
      SplitPattern_2_2_9_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_2_9_Dispatcher{{weight[0]}, {weight[1], delta[0]}, {weight[2], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {1.105332, 1.143884};
    static constexpr std::array<float, 3> weight = {.213636, .126689, .069902};
};


/**************************************************************************/
using SplitPattern_2_2_16_Dispatcher
    = PatternDispatcher<PurplePattern<DimConst<2>>, BrownPattern<DimConst<2>>,
                        PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<16>>
    : public ASplitter<DimConst<2>, InterpConst<2>, RefinedParticlesConst<16>>,
      SplitPattern_2_2_16_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_2_16_Dispatcher({weight[0], delta[0]}, {weight[1], delta},
                                         {weight[2], delta[1]})
    {
    }

    static constexpr std::array<float, 2> delta  = {.5, 1.5};
    static constexpr std::array<float, 3> weight = {0.140625, 0.046875, 0.015625};
};


/**************************************************************************/
using SplitPattern_2_3_4_Dispatcher = PatternDispatcher<PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<4>>
    : public ASplitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<4>>,
      SplitPattern_2_3_4_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_3_4_Dispatcher{{weight[0], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {0.776459};
    static constexpr std::array<float, 1> weight = {0.25};
};


/**************************************************************************/
using SplitPattern_2_3_5_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PinkPattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<5>>
    : public ASplitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<5>>,
      SplitPattern_2_3_5_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_3_5_Dispatcher{{weight[0]}, {weight[1], delta[0]}}
    {
    }

    static constexpr std::array<float, 1> delta  = {1.376953};
    static constexpr std::array<float, 2> weight = {0.242666, 0.189333};
};


/**************************************************************************/
using SplitPattern_2_3_8_Dispatcher
    = PatternDispatcher<PinkPattern<DimConst<2>>, PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<8>>
    : public ASplitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<8>>,
      SplitPattern_2_3_8_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_3_8_Dispatcher{{weight[0], delta[0]}, {weight[1], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {0.942365, 1.423324};
    static constexpr std::array<float, 2> weight = {0.179318, 0.070682};
};


/**************************************************************************/
using SplitPattern_2_3_9_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PinkPattern<DimConst<2>>,
                        PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<9>>
    : public ASplitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<9>>,
      SplitPattern_2_3_9_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_3_9_Dispatcher{{weight[0]}, {weight[1], delta[0]}, {weight[2], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta  = {1.267689, 1.315944};
    static constexpr std::array<float, 3> weight = {.218605, .126871, .068477};
};


/**************************************************************************/
using SplitPattern_2_3_25_Dispatcher
    = PatternDispatcher<BlackPattern<DimConst<2>>, PinkPattern<DimConst<2>>,
                        PurplePattern<DimConst<2>>, PinkPattern<DimConst<2>>,
                        BrownPattern<DimConst<2>>, PurplePattern<DimConst<2>>>;

template<>
struct Splitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<25>>
    : public ASplitter<DimConst<2>, InterpConst<3>, RefinedParticlesConst<25>>,
      SplitPattern_2_3_25_Dispatcher
{
    constexpr Splitter()
        : SplitPattern_2_3_25_Dispatcher{{weight[0]},           {weight[1], delta[0]},
                                         {weight[2], delta[0]}, {weight[3], delta[1]},
                                         {weight[4], delta},    {weight[5], delta[1]}}
    {
    }

    static constexpr std::array<float, 2> delta = {1, 2};
    static constexpr std::array<float, 6> weight
        = {.140625, .09375, .0625, .0234375, .015625, .00390625};
};


/**************************************************************************/


} // namespace PHARE::amr


#endif /*PHARE_SPLIT_2D_H*/
