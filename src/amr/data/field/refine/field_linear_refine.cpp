#include "field_linear_refine.h"

using namespace PHARE;

UniformIntervalPartitionWeight::UniformIntervalPartitionWeight(QtyCentering centering, int ratio,
                                                               std::size_t nbrPoints)
{
    assert(nbrPoints > 1);

    distances_.resize(nbrPoints);

    bool evenRatio = ratio % 2 == 0;


    // when we are primal we have the coarse centering
    // that lie on top of a fine centered data
    if (centering == QtyCentering::primal)
    {
        for (std::size_t i = 0; i < distances_.size(); ++i)
        {
            distances_[i] = static_cast<double>(i) / ratio;
        }
    }
    else
    {
        if (evenRatio)
        {
            int const halfRatio = ratio / 2;
            for (std::size_t i = 0; i < distances_.size(); ++i)
            {
                int const j = (halfRatio + i) % ratio;

                distances_[j] = 1 / (2. * ratio) + static_cast<double>(i) / ratio;
            }
        }
        // when not evenRatio, we have  a coarse centering  that lie on top of a fine centered data
        else
        {
            int const halfRatio = ratio / 2;
            for (std::size_t i = 0; i < distances_.size(); ++i)
            {
                int const j = (halfRatio + i) % ratio;

                distances_[j] = static_cast<double>(i) / ratio;
            }
        }
    }
}
