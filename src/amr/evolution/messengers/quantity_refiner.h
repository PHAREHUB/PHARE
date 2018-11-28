#ifndef PHARE_QUANTITY_REFINER_H
#define PHARE_QUANTITY_REFINER_H


#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>

#include <map>
#include <memory>
#include <optional>



/**
 * @brief The Refiner struct encapsulate the algorithm and its associated
 * schedules
 *
 * We have several of those object, one per level, that is used to retrieve which schedule to
 * use for a given messenger communication, and which algorithm to use to re-create schedules
 * when initializing a level
 */
struct QuantityRefiner
{
    QuantityRefiner()
        : algo{std::make_unique<SAMRAI::xfer::RefineAlgorithm>()}
    {
    }

    std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> algo; // this part is set in setup

    // this part is created in initializeLevelData()
    std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> schedules;

    std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>> findSchedule(int levelNumber)
    {
        if (auto mapIter = schedules.find(levelNumber); mapIter != std::end(schedules))
        {
            return mapIter->second;
        }
        else
        {
            return std::nullopt;
        }
    }

    void add(std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule, int levelNumber)
    {
        schedules[levelNumber] = std::move(schedule);
    }
};


struct RefinerPool
{
    std::map<std::string, QuantityRefiner> qtyRefiners;

    std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
    findSchedule(std::string const& name, int levelNumber)
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


    void add(QuantityRefiner&& qtyRefiner, std::string const& qtyName)
    {
        if (qtyRefiners.find(qtyName) == std::end(qtyRefiners))
            qtyRefiners[qtyName] = std::move(qtyRefiner);
        else
            throw std::runtime_error(qtyName + " already in map");
    }
};




#endif
