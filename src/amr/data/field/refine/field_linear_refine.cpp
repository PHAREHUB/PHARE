#include "field_linear_refine.h"



namespace PHARE
{
UniformIntervalPartitionWeight::UniformIntervalPartitionWeight(QtyCentering centering,
                                                               std::size_t ratio,
                                                               std::size_t nbrPoints)

{
    assert(nbrPoints > 1);
    distances_.resize(nbrPoints);
    bool isEvenRatio   = ratio % 2 == 0;
    auto smallCellSize = 1. / ratio;

    std::iota(std::begin(distances_), std::end(distances_), 0);

    // when we are primal we have the coarse centering
    // that lie on top of a fine centered data
    if (centering == QtyCentering::primal)
    {
        std::transform(std::begin(distances_), std::end(distances_), std::begin(distances_),
                       [ratio](auto const& v) { return static_cast<double>(v) / ratio; });
    }
    else
    {
        if (isEvenRatio)
        {
            auto middle = std::begin(distances_) + distances_.size() / 2;
            std::transform(std::begin(distances_), std::end(distances_), std::begin(distances_),
                           [smallCellSize](auto const& v) {
                               return (0.5 + static_cast<double>(v)) * smallCellSize;
                           });
            std::rotate(std::begin(distances_), middle, std::end(distances_));
        }

        else
        {
            auto middle = std::begin(distances_) + distances_.size() / 2 + 1;
            std::transform(
                std::begin(distances_), std::end(distances_), std::begin(distances_),
                [smallCellSize](auto const& v) { return static_cast<double>(v) * smallCellSize; });

            std::rotate(std::begin(distances_), middle, std::end(distances_));
        }
    }
}

} // namespace PHARE
