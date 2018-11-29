

#include "quantity_refiner.h"

namespace PHARE
{
std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
RefinerPool::findSchedule_(std::string const& name, int levelNumber)
{
    if (auto mapIter = qtyRefiners.find(name); mapIter != std::end(qtyRefiners))
    {
        return mapIter->second.findSchedule(levelNumber);
    }
    else
    {
        return std::nullopt;
    }
}




void RefinerPool::add(QuantityRefiner&& qtyRefiner, std::string const& qtyName)
{
    if (qtyRefiners.find(qtyName) == std::end(qtyRefiners))
        qtyRefiners[qtyName] = std::move(qtyRefiner);
    else
        throw std::runtime_error(qtyName + " already in map");
}



void RefinerPool::createGhostSchedules(
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
    std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
{
    for (auto& [key, refiner] : qtyRefiners)
    {
        auto& algo = refiner.algo;
        auto schedule
            = algo->createSchedule(level, level->getNextCoarserHierarchyLevelNumber(), hierarchy);
        refiner.add(schedule, level->getLevelNumber());
    }
}




void RefinerPool::createInitSchedules(
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
    std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
{
    for (auto& [key, refiner] : qtyRefiners)
    {
        auto& algo       = refiner.algo;
        auto levelNumber = level->getLevelNumber();

        // note that here we must take that createsSchedule() overload and put nullptr as src
        // since we want to take from coarser level everywhere.
        // using the createSchedule overload that takes level, next_coarser_level only
        // would result in interior ghost nodes to be filled with interior of neighbor patches
        // but there is nothing there.
        refiner.add(algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy), levelNumber);
    }
}




void RefinerPool::initialize(int levelNumber, double initDataTime) const
{
    for (auto& [key, refiner] : qtyRefiners)

    {
        if (refiner.algo == nullptr)
        {
            throw std::runtime_error("Algorithm is nullptr");
        }

        auto schedule = refiner.schedules.find(levelNumber);
        if (schedule != std::end(refiner.schedules))
        {
            schedule->second->fillData(initDataTime);
        }
        else
        {
            throw std::runtime_error("Error - schedule cannot be found for this level");
        }
    }
}


} // namespace PHARE
