#ifndef PHARE_AMR_SOLVERS_PARTICLE_LEVEL_BOXINGS_HPP
#define PHARE_AMR_SOLVERS_PARTICLE_LEVEL_BOXINGS_HPP

#include "core/numerics/ion_updater/ion_updater.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>


namespace PHARE::solver
{
/**
 * @brief ParticleLevelBoxings keeps, per level and per patch, the boxes the particle updater
 * selects particles with (domain, owned cells, non level ghost box). They are built the first
 * time a level is advanced and rebuilt when the level object changes, even without a regrid
 * notification.
 */
template<typename IonUpdater_t, typename GridLayout, typename AMR_Types>
class ParticleLevelBoxings
{
    using hierarchy_t = AMR_Types::hierarchy_t;
    using level_t     = AMR_Types::level_t;
    using level_ptr   = AMR_Types::level_ptr;

public:
    using Boxing_t = core::UpdaterSelectionBoxing<IonUpdater_t, GridLayout>;

    struct LevelBoxing
    {
        std::weak_ptr<level_t> level; // identity of the level the boxes were built for
        std::unordered_map<std::string /*patchid*/, Boxing_t> patches;
    };

    auto& setup_level(hierarchy_t const& hierarchy, int const levelNumber)
    {
        auto level = hierarchy.getPatchLevel(levelNumber);
        // level may have been removed and re-created without a regrid notification
        auto const it            = boxing_.find(levelNumber);
        bool const needsCreating = it == boxing_.end() or it->second.level.lock() != level;
        if (needsCreating)
            make_boxes_(hierarchy, level);
        return *level;
    }

    auto& operator[](int const levelNumber) { return boxing_[levelNumber]; }

    void clear() { boxing_.clear(); }

private:
    void make_boxes_(hierarchy_t const& hierarchy, level_ptr const& level)
    {
        auto& levelBoxing = boxing_[level->getLevelNumber()] = LevelBoxing{level, {}};

        auto const neighbors = amr::makeSameLevelNeighbors(hierarchy, level->getLevelNumber());

        for (auto const& patch : *level)
            if (auto [it, suc] = levelBoxing.patches.try_emplace(
                    amr::to_string(patch->getGlobalId()),
                    Boxing_t{amr::layoutFromPatch<GridLayout>(*patch),
                             amr::makeNonLevelGhostBoxFor<GridLayout>(*patch, neighbors),
                             amr::makeNonOwnedBoxesFor<GridLayout::dimension>(*patch, neighbors)});
                !suc)
                throw std::runtime_error("boxing map insertion failure");
    }

    std::unordered_map<int /*level*/, LevelBoxing> boxing_;
};

} // namespace PHARE::solver

#endif
