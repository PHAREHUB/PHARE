#ifndef PHARE_LINEAR_WEIGHTER_H
#define PHARE_LINEAR_WEIGHTER_H


#include "data/field/field.h"
#include "data/grid/gridlayoutdefs.h"
#include "utilities/constants.h"
#include "utilities/point/point.h"

#include <SAMRAI/hier/Box.h>

#include <array>
#include <vector>


namespace PHARE
{
namespace amr_interface
{
    /** @brief  This class calculates the distances of each fine index within a coarse cell
     * from the left-most coarse index of the same kind (dual/primal)
     */
    class LinearWeighter
    {
    public:
        using FineIndexWeight  = std::array<double, 2>;
        using FineIndexWeights = std::vector<FineIndexWeight>;


        LinearWeighter(core::QtyCentering centering, std::size_t ratio);

        std::vector<double> const& getUniformDistances() const { return distFromLeftNode_; }
        FineIndexWeights const& weights() const { return weights_; }

    private:
        std::vector<double> distFromLeftNode_;
        FineIndexWeights weights_;
    };


    template<typename T, std::size_t... Is>
    std::array<LinearWeighter, sizeof...(Is)>
    make_weighters(const std::array<T, sizeof...(Is)>& values, SAMRAI::hier::IntVector ratio,
                   std::index_sequence<Is...>)
    {
        return {{(LinearWeighter{values[Is], static_cast<std::size_t>(ratio[Is])})...}};
    }
} // namespace amr_interface
} // namespace PHARE


#endif
