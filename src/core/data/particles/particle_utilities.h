#ifndef PHARE_PARTICLE_UTILITIES
#define PHARE_PARTICLE_UTILITIES

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/particles/particle.h"
#include "core/utilities/point/point.h"

#include <array>


namespace PHARE::core
{
template<typename GridLayout>
/**
 * @brief positionAsPoint returns a point holding the physical position of the macroparticle.
 * The function assumes the iCell of the particle is in AMR index space.
 */
auto positionAsPoint(Particle<GridLayout::dimension> const& particle, GridLayout const& layout)
{
    Point<double, GridLayout::dimension> position;
    auto origin       = layout.origin();
    auto startIndexes = layout.physicalStartIndex(QtyCentering::primal);
    auto meshSize     = layout.meshSize();
    auto iCell        = layout.AMRToLocal(Point<int, GridLayout::dimension>{particle.iCell});

    for (auto iDim = 0u; iDim < GridLayout::dimension; ++iDim)
    {
        position[iDim] = origin[iDim];
        position[iDim]
            += (iCell[iDim] - startIndexes[iDim] + particle.delta[iDim]) * meshSize[iDim];
    }
    return position;
}


// this function is for debugging purposes
// it looks for the given particle range, whether two particles
// are found to have the same "icell+delta" (in X direction)
// this is not a problem in theory but makes the test
// test_overlapped_particledatas_have_identical_particles
// fail
template<typename ParticleRange>
void checkDeltas(ParticleRange const& prange)
{
    std::vector<double> deltas;

    for (auto const& part : prange)
    {
        deltas.push_back(part.iCell[0] + part.delta[0]);
    }
    std::sort(std::begin(deltas), std::end(deltas));
    auto p = std::adjacent_find(std::begin(deltas), std::end(deltas));
    if (p != std::end(deltas))
    {
        double delta = *p;
        std::cout << "Oops found duplicates \n";
        auto part = std::find_if(std::begin(prange), std::end(prange), [delta](auto const& par) {
            return par.iCell[0] + par.delta[0] == delta;
        });

        auto part2 = std::find_if(part + 1, std::end(prange), [delta](auto const& par) {
            return par.iCell[0] + par.delta[0] == delta;
        });

        if (part == std::end(prange))
            std::cout << "part at the end of prange\n";
        else if (part2 == std::end(prange))
            std::cout << "part2 at the end of prange\n";
        else
        {
            std::cout << std::setprecision(12) << "part1 delta = " << part->delta[0]
                      << " , part2 delta = " << part2->delta[0] << "\n";
            std::cout << "part1 cell = " << part->iCell[0] << " , part2 Cell = " << part2->iCell[0]
                      << "\n";
            std::cout << "distance : " << std::distance(part, part2) << "\n";
            std::cout << "size : " << prange.size() << "\n";
            std::cout << "p1.vx = " << part->v[0] << " p2.vx = " << part2->v[0] << "\n";
            std::cout << "p1.vy = " << part->v[1] << " p2.vy = " << part2->v[1] << "\n";
            std::cout << "p1.vz = " << part->v[2] << " p2.vz = " << part2->v[2] << "\n";
            std::cout << &(*part) << " " << &(*part2) << "\n";
            throw std::runtime_error("Oops duplicates duplicates throooow");
        }
    }
    else
        std::cout << "no duplicate found\n";
}



} // namespace PHARE::core



#endif
