

#include "quantity_refiner.h"

namespace PHARE
{
std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
RefinerPool::findSchedule_(std::string const& name, int levelNumber)
{
    if (auto mapIter = qtyRefiners_.find(name); mapIter != std::end(qtyRefiners_))
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
    if (qtyRefiners_.find(qtyName) == std::end(qtyRefiners_))
        qtyRefiners_[qtyName] = std::move(qtyRefiner);
    else
        throw std::runtime_error(qtyName + " already in map");
}



void RefinerPool::createGhostSchedules(
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
    std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
{
    for (auto& [key, refiner] : qtyRefiners_)
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
    for (auto& [key, refiner] : qtyRefiners_)
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
    for (auto& [key, refiner] : qtyRefiners_)
    {
        if (refiner.algo == nullptr)
        {
            throw std::runtime_error("Algorithm is nullptr");
        }

        auto schedule = refiner.findSchedule(levelNumber);
        if (schedule)
        {
            (*schedule)->fillData(initDataTime);
        }
        else
        {
            throw std::runtime_error("Error - schedule cannot be found for this level");
        }
    }
}




void RefinerPool::regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                         const int levelNumber,
                         std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                         double const initDataTime)
{
    for (auto& [key, refiner] : qtyRefiners_)
    {
        auto& algo = refiner.algo;

        // here 'nullptr' is for 'oldlevel' which is always nullptr in this function
        // the regriding schedule for which oldptr is not nullptr is handled in another
        // function
        auto const& level = hierarchy->getPatchLevel(levelNumber);

        auto schedule = algo->createSchedule(
            level, oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy);

        schedule->fillData(initDataTime);
    }
}




} // namespace PHARE
