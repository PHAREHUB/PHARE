#ifndef PHARE_PARTICLE_UTILITIES
#define PHARE_PARTICLE_UTILITIES

#include "data/grid/gridlayoutdefs.h"
#include "data/particles/particle.h"
#include "utilities/point/point.h"

#include <array>


namespace PHARE
{
template<typename GridLayout>
auto positionAsPoint(Particle<GridLayout::dimension> const& particle, GridLayout const& layout)
{
    Point<double, GridLayout::dimension> position;
    auto origin       = layout.origin();
    auto startIndexes = layout.physicalStartIndex(QtyCentering::primal);
    auto meshSize     = layout.meshSize();


    for (auto iDim = 0u; iDim < GridLayout::dimension; ++iDim)
    {
        position[iDim] = origin[iDim];
        position[iDim]
            += (particle.iCell[iDim] - startIndexes[iDim] + particle.delta[iDim]) * meshSize[iDim];
    }
    return position;
}




} // namespace PHARE



#endif
