#include "field_coarsen_index_weight.hpp"

namespace PHARE::amr
{
[[nodiscard]] double CoarsenWeighter::findX_(std::size_t nbrPoints) const
{
    double x = 0.;

    if (nbrPoints % 2 != 0)
    {
        x = 1.;
        for (std::size_t i = 1; i <= (nbrPoints - 1) / 2; ++i)
        {
            x += 2 * 1. / (i + 1);
        }
    }
    else
    {
        for (std::size_t i = 1; i <= nbrPoints / 2; ++i)
        {
            x += 2 * 1. / i;
        }
    }

    return x;
}




void CoarsenWeighter::computeWeights_(std::size_t nbrPoints)
{
    weights_.resize(nbrPoints);

    auto x = findX_(nbrPoints);


    if (nbrPoints % 2 != 0)
    {
        auto const halfIndex = (nbrPoints - 1) / 2;

        auto const halfNumberOfPointsLeft
            = nbrPoints / 2; // half of the points needed besides the one on the middle

        weights_[halfIndex] = 1. / x;

        for (std::size_t i = 1; i <= halfNumberOfPointsLeft; ++i)
        {
            double factor           = static_cast<double>(i + 1);
            weights_[halfIndex - i] = 1. / (factor * x);
            weights_[halfIndex + i] = 1. / (factor * x);
        }
    }




    else
    {
        auto const halfIndexRight = nbrPoints / 2;
        auto const halfIndexLeft  = halfIndexRight - 1;

        weights_[halfIndexRight] = 1. / x;
        weights_[halfIndexLeft]  = 1. / x;

        auto const halfNumberOfPointsLeft = (nbrPoints / 2) - 1;

        for (std::size_t i = 1; i <= halfNumberOfPointsLeft; ++i)
        {
            double factor = static_cast<double>(i + 1);

            weights_[halfIndexRight + i] = 1. / (factor * x);
            weights_[halfIndexLeft - i]  = 1. / (factor * x);
        }
    }
}


} // namespace PHARE::amr
