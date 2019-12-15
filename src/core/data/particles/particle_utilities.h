#ifndef PHARE_PARTICLE_UTILITIES
#define PHARE_PARTICLE_UTILITIES

#include "core/data/grid/gridlayoutdefs.h"
#include "core/utilities/point/point.h"
#include "core/data/particles/particle.h"

#include <array>


namespace PHARE
{
namespace core
{
    /**
     * @brief positionAsPoint returns a point holding the physical position of the macroparticle.
     * The function assumes the iCell of the particle is in AMR index space.
     */
    template<typename Particle, typename GridLayout>
    auto positionAsPoint(Particle const& particle, GridLayout const& layout)
    {
        Point<double, GridLayout::dimension> position;
        auto origin       = layout.origin();
        auto startIndexes = layout.physicalStartIndex(QtyCentering::primal);
        auto meshSize     = layout.meshSize();
        auto iCell        = layout.AMRToLocal(Point{particle.iCell});

        for (auto iDim = 0u; iDim < GridLayout::dimension; ++iDim)
        {
            position[iDim] = origin[iDim];
            position[iDim]
                += (iCell[iDim] - startIndexes[iDim] + particle.delta[iDim]) * meshSize[iDim];
        }
        return position;
    }
} // namespace core




} // namespace PHARE



#endif
