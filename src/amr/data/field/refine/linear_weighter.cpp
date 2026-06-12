#include "field_linear_refine.hpp"

#include <algorithm>

namespace PHARE::amr
{
LinearWeighter::LinearWeighter(core::QtyCentering centering, std::size_t ratio)

{
    auto nbrPoints = ratio;
    assert(nbrPoints > 1);
    distFromLeftNode_.resize(nbrPoints);
    bool isEvenRatio   = ratio % 2 == 0;
    auto smallCellSize = 1. / ratio;

    std::iota(std::begin(distFromLeftNode_), std::end(distFromLeftNode_), 0);

    // when we are primal we have the coarse centering
    // that lie on top of a fine centered data
    if (centering == core::QtyCentering::primal)
    {
        std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                       std::begin(distFromLeftNode_),
                       [ratio](auto const& v) { return static_cast<double>(v) / ratio; });
    }
    else
    {
        if (isEvenRatio)
        {
            auto middle = std::begin(distFromLeftNode_) + distFromLeftNode_.size() / 2;
            std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                           std::begin(distFromLeftNode_), [smallCellSize](auto const& v) {
                               return (0.5 + static_cast<double>(v)) * smallCellSize;
                           });
            std::rotate(std::begin(distFromLeftNode_), middle, std::end(distFromLeftNode_));
        }

        else
        {
            auto middle = std::begin(distFromLeftNode_) + distFromLeftNode_.size() / 2 + 1;
            std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                           std::begin(distFromLeftNode_), [smallCellSize](auto const& v) {
                               return static_cast<double>(v) * smallCellSize;
                           });

            std::rotate(std::begin(distFromLeftNode_), middle, std::end(distFromLeftNode_));
        }
    }


    std::transform(std::begin(distFromLeftNode_), std::end(distFromLeftNode_),
                   std::back_inserter(weights_),
                   [](auto const& d) { return std::array<double, 2>{{1. - d, d}}; });
}

} // namespace PHARE::amr
