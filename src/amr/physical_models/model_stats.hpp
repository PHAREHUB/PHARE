#ifndef PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP
#define PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP


#include "core/utilities/types.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/resources_manager/amr_utils.hpp"


#include <ostream>


namespace PHARE::solver
{
template<typename GridLayoutT, typename PhysicalQuantity>
struct LevelStats
{
    auto static constexpr primal_qty = PhysicalQuantity::all_primal_field;

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
            auto const box       = layout.AMRBoxFor(primal_qty);
            auto const ghost_box = layout.AMRGhostBoxFor(primal_qty);

            auto const cells  = box.size();
            auto const ghosts = ghost_box.size() - cells;
            auto const patch_ghost_cells
                = sum_from(amr::patchGhostBoxOverlaps<GridLayoutT>(*patch, neighbors, ghost_width),
                           [](auto box) {
                               box.upper += 1;
                               return box.size();
                           });

            stats.cells += cells;
            stats.patch_ghosts += patch_ghost_cells;
            stats.level_ghosts += ghosts - patch_ghost_cells;
            ++stats.patches;
        }

        return stats;
    }

    auto collect() const
    {
        using namespace PHARE::core;
        LevelStats stats;
        stats.cells        = mpi::sum_on_rank_0(cells);
        stats.patch_ghosts = mpi::sum_on_rank_0(patch_ghosts);
        stats.level_ghosts = mpi::sum_on_rank_0(level_ghosts);
        stats.patches      = mpi::sum_on_rank_0(patches);
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
