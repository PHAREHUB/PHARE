#ifndef PHARE_COARSEN_WEIGHTER_H
#define PHARE_COARSEN_WEIGHTER_H

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>


namespace PHARE
{
namespace amr
{
    /**
     * @brief The CoarsenWeighter class computes the weights to use for each nbrPoints fine nodes to
     * get the value on the associated coarse node
     */
    class CoarsenWeighter
    {
    public:
        explicit CoarsenWeighter(std::size_t nbrPoints)
        {
            assert(nbrPoints > 1); // we want to have at least 2 points for coarsening operations
            computeWeights_(nbrPoints);
        }

        std::vector<double> const& weights() const { return weights_; }

    private:
        std::vector<double> weights_;

        double findX_(std::size_t nbrPoints) const;
        void computeWeights_(std::size_t nbrPoints);
    };



    template<typename T, std::size_t... Is>
    std::array<CoarsenWeighter, sizeof...(Is)>
    make_weighters(const std::array<T, sizeof...(Is)>& nbrPoints, std::index_sequence<Is...>)
    {
        return {{(CoarsenWeighter{nbrPoints[Is]})...}};
    }
} // namespace amr

} // namespace PHARE

#endif
