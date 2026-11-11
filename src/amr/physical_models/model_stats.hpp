#ifndef PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP
#define PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP


#include "mpi/mpi_utils.hpp"

#include "amr/resources_manager/amr_utils.hpp"


#include <map>
#include <set>
#include <string>
#include <cstdint>
#include <variant>
#include <stdexcept>

#include <SAMRAI/hier/Connector.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/BoxContainer.h>


#if !defined(PHARE_SUMMARY_MODE)
#define PHARE_SUMMARY_MODE 0 // 0 light, 1 heavy: adds SAMRAI load balance reports
#endif

namespace PHARE
{
// light: summary() only reports, no side effects outside the call
// heavy: also enables diagnostics that run on their own (e.g. on each load balance)
constexpr static std::uint8_t SUMMARY_MODE = PHARE_SUMMARY_MODE;
} // namespace PHARE


namespace PHARE::solver
{
// one level of Simulator::summary(), python gets a dict of whichever type each value holds
using LevelSummary = std::map<std::string, std::variant<std::size_t, double, std::string>>;


// max per rank over average per rank, 1 is perfectly balanced
inline double imbalance(std::size_t const max, std::size_t const total)
{
    return total ? static_cast<double>(max) * mpi::size() / total : 1;
}


template<typename GridLayoutT, typename PhysicalQuantity>
struct LevelStats
{
    std::size_t patches = 0, cells = 0, patch_ghosts = 0, level_ghosts = 0;
    std::size_t remote_ghosts = 0; // the patch ghosts owned by another rank
    std::size_t remote_ranks  = 0; // other ranks owning patch ghosts, max over ranks once collected
    double remote_ranks_avg   = 0; // average over ranks, once collected
    double cell_imbalance     = 1; // max cells per rank over average, 1 is perfectly balanced

    // ghost cells within ghost_width of the patch box: level ghosts are those not covered by
    // the patch or any same level neighbour, patch ghosts are the rest, split by neighbour owner
    static auto GET(auto& lvl, int const ghost_width)
    {
        auto const& connector = sameLevelConnector_(lvl, ghost_width);
        auto const rank       = mpi::rank();

        LevelStats stats;
        std::set<int> remote_ranks;
        for (auto const& patch : lvl)
        {
            auto const& domain = patch->getBox();
            auto ghost_box     = domain;
            ghost_box.grow(SAMRAI::hier::IntVector{domain.getDim(), ghost_width});

            // includes the patch itself and periodic images, already shifted
            SAMRAI::hier::BoxContainer neighbors;
            if (connector.hasNeighborSet(domain.getBoxId()))
                connector.getNeighborBoxes(domain.getBoxId(), neighbors);

            // exact set differences: overlapping neighbours are not counted twice
            SAMRAI::hier::BoxContainer level_ghost_boxes{ghost_box};
            level_ghost_boxes.removeIntersections(domain);
            for (auto const& neighbor : neighbors)
                if (neighbor.getOwnerRank() == rank)
                    level_ghost_boxes.removeIntersections(neighbor);
            auto const not_local = level_ghost_boxes.getTotalSizeOfBoxes();
            for (auto const& neighbor : neighbors)
                if (neighbor.getOwnerRank() != rank)
                {
                    level_ghost_boxes.removeIntersections(neighbor);
                    remote_ranks.insert(neighbor.getOwnerRank());
                }

            auto const level_ghosts = level_ghost_boxes.getTotalSizeOfBoxes();
            stats.level_ghosts += level_ghosts;
            stats.remote_ghosts += not_local - level_ghosts;
            stats.patch_ghosts += ghost_box.size() - domain.size() - level_ghosts;
        }
        stats.remote_ranks = remote_ranks.size();

        return stats;
    }

    // patches, cells and max cells per rank from SAMRAI's cached global reduction, collective on
    // first call
    auto collect(auto const& lvl) const
    {
        LevelStats stats;
        stats.patches        = lvl.getGlobalNumberOfPatches();
        stats.cells          = lvl.getGlobalNumberOfCells();
        stats.cell_imbalance = imbalance(lvl.getBoxLevel()->getMaxNumberOfCells(), stats.cells);
        stats.patch_ghosts   = mpi::sum_on_rank_0(patch_ghosts);
        stats.level_ghosts   = mpi::sum_on_rank_0(level_ghosts);
        stats.remote_ghosts  = mpi::sum_on_rank_0(remote_ghosts);
        stats.remote_ranks   = mpi::max(remote_ranks);
        stats.remote_ranks_avg
            = static_cast<double>(mpi::sum_on_rank_0(remote_ranks)) / mpi::size();
        return stats;
    }

    // collected values, models add theirs, e.g. "particles"
    LevelSummary summary(std::size_t const level, std::string const& model) const
    {
        return {
            {"level", level},
            {"model", model},
            {"patches", patches},
            {"cells", cells},
            {"patch_ghosts", patch_ghosts},
            {"level_ghosts", level_ghosts},
            {"remote_ghosts", remote_ghosts},
            {"remote_ranks_max", remote_ranks},
            {"remote_ranks_avg", remote_ranks_avg},
            {"cell_imbalance", cell_imbalance},
        };
    }

private:
    // the existing same level connector at least ghost_width wide, never creates one: light
    // summaries have no side effects. The patch ghost refine schedules need one at least as wide
    // as the field ghost width.
    static auto& sameLevelConnector_(auto const& lvl, int const ghost_width)
    {
        SAMRAI::hier::IntVector const width{lvl.getDim(), ghost_width};
        if (!lvl.hasConnector(lvl, width))
            throw std::runtime_error("LevelStats: no same level connector of width "
                                     + std::to_string(ghost_width));
        return lvl.findConnector(lvl, width, SAMRAI::hier::CONNECTOR_ERROR,
                                 /*exact_width_only*/ false);
    }
};

} // namespace PHARE::solver

#endif /* PHARE_AMR_PHYSICAL_MODELS_MODEL_STATS_HPP */
