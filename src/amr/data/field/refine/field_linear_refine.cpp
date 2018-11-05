#include "field_linear_refine.h"



namespace PHARE
{
std::size_t nbrRefinePoints(QtyCentering centering, bool isEvenRatio, int ratio)
{
    std::size_t nbrPoints;
    if (centering == QtyCentering::primal || !isEvenRatio)
    {
        nbrPoints = static_cast<std::size_t>(ratio);
    }
    else
    {
        nbrPoints = static_cast<std::size_t>(ratio);
    }
    return nbrPoints;
}


UniformIntervalPartitionWeight::UniformIntervalPartitionWeight(QtyCentering centering, int ratio,
                                                               std::size_t nbrPoints)
{
    assert(nbrPoints > 1);

    distances_.resize(nbrPoints);

    bool isEvenRatio   = ratio % 2 == 0;
    auto smallCellSize = 1. / ratio;

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
        if (isEvenRatio)
        {
            int const halfRatio = ratio / 2;
            for (std::size_t i = 0; i < distances_.size(); ++i)
            {
                int const j = (halfRatio + i) % ratio;

                distances_[j] = (1 / 2. + static_cast<double>(i)) * smallCellSize;
            }
        }
        // when not evenRatio, we have  a coarse centering  that lie on top of a fine centered data
        else
        {
            int const halfRatio = ratio / 2;
            for (std::size_t i = 0; i < distances_.size() - 1; ++i)
            {
                int const j = (halfRatio + i) % ratio;

                distances_[j] = static_cast<double>(i) / ratio;
            }
        }
    }
}

} // namespace PHARE
