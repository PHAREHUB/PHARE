#ifndef PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP
#define PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP


#include "core/utilities/box/box_spans.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/resources_manager/amr_utils.hpp"


#include <ostream>


namespace PHARE::solver
{
template<typename GridLayoutT, typename PhysicalQuantity>
struct LevelStats
{
    std::size_t patches = 0, cells = 0, patch_ghosts = 0, level_ghosts = 0;

    // patch ghosts are counted ghost_width cells out from the patch box
    static auto GET(auto& lvl, auto& hierarchy, int const ghost_width)
    {
        using namespace PHARE::core;

        auto const neighbors = amr::makeSameLevelNeighbors(hierarchy, lvl.getLevelNumber());

        LevelStats stats;
        for (auto const& patch : lvl)
        {
            auto const layout    = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto const box       = layout.AMRBox();
            auto const ghost_box = grow(box, ghost_width);

            auto const ghosts = ghost_box.size() - box.size();

            // neighbours' cells in the ghost box, each ghost cell counted once: overlaps can
            // overlap each other
            auto const neighbour_cells
                = amr::patchGhostBoxOverlaps<GridLayoutT>(*patch, neighbors, ghost_width);
            std::size_t patch_ghost_cells = 0;
            for_each_point(ghost_box, [&](auto const& cell) {
                if (!isIn(cell, box) and isIn(cell, neighbour_cells))
                    ++patch_ghost_cells;
            });

            stats.patch_ghosts += patch_ghost_cells;
            stats.level_ghosts += ghosts - patch_ghost_cells;
        }

        return stats;
    }

    // patches and cells from SAMRAI's cached global reduction, collective on first call
    auto collect(auto const& lvl) const
    {
        using namespace PHARE::core;
        LevelStats stats;
        stats.patches      = lvl.getGlobalNumberOfPatches();
        stats.cells        = lvl.getGlobalNumberOfCells();
        stats.patch_ghosts = mpi::sum_on_rank_0(patch_ghosts);
        stats.level_ghosts = mpi::sum_on_rank_0(level_ghosts);
        return stats;
    }
};

template<typename GridLayoutT, typename PhysicalQuantity>
inline std::ostream& operator<<(std::ostream& out,
                                LevelStats<GridLayoutT, PhysicalQuantity> const& stats)
{
    out << stats.patches << ",";
    out << stats.cells << ",";
    out << stats.patch_ghosts << ",";
    out << stats.level_ghosts;

    return out;
}


} // namespace PHARE::solver

#endif /* PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP */
