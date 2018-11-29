#ifndef PHARE_QUANTITY_REFINER_H
#define PHARE_QUANTITY_REFINER_H


#include "evolution/messengers/hybrid_messenger_info.h"

#include "data/field/refine/field_refine_operator.h"
#include "data/field/time_interpolate/field_linear_time_interpolate.h"
#include "data/particles/refine/particles_data_split.h"
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>

#include <map>
#include <memory>
#include <optional>

namespace PHARE
{
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




template<typename ResourcesManager>
QuantityRefiner makeGhostRefiner(VecFieldDescriptor const& ghost, VecFieldDescriptor const& model,
                                 VecFieldDescriptor const& oldModel, ResourcesManager const& rm,
                                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                                 std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp)
{
    QuantityRefiner refiner;
    refiner.algo = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();

    auto registerRefine
        = [&rm, &refiner, &refineOp, &timeOp](std::string const& model, std::string const& ghost,
                                              std::string const& oldModel) {
              auto src_id  = rm->getID(model);
              auto dest_id = rm->getID(ghost);
              auto old_id  = rm->getID(oldModel);

              if (src_id && dest_id && old_id)
              {
                  // dest, src, old, new, scratch
                  refiner.algo->registerRefine(*dest_id, // dest
                                               *src_id,  // source at same time
                                               *old_id,  // source at past time (for time interp)
                                               *src_id,  // source at future time (for time interp)
                                               *dest_id, // scratch
                                               refineOp, timeOp);
              }
          };

    // register refine operators for each component of the vecfield
    registerRefine(ghost.xName, model.xName, oldModel.xName);
    registerRefine(ghost.yName, model.yName, oldModel.yName);
    registerRefine(ghost.zName, model.zName, oldModel.zName);

    return refiner;
}



template<typename ResourcesManager>
QuantityRefiner makeInitRefiner(VecFieldDescriptor const& name, ResourcesManager const& rm,
                                std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
{
    QuantityRefiner refiner;
    refiner.algo = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();

    auto registerRefine = [&refiner, &rm, &refineOp](std::string name) //
    {
        auto id = rm->getID(name);
        if (id)
        {
            refiner.algo->registerRefine(*id, *id, *id, refineOp);
        }
    };

    registerRefine(name.xName);
    registerRefine(name.yName);
    registerRefine(name.zName);

    return refiner;
}




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



    void createGhostSchedules(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
    {
        for (auto& [key, refiner] : qtyRefiners)
        {
            auto& algo    = refiner.algo;
            auto schedule = algo->createSchedule(level, level->getNextCoarserHierarchyLevelNumber(),
                                                 hierarchy);
            refiner.add(schedule, level->getLevelNumber());
        }
    }


    template<typename VecFieldT>
    void fillVecFieldGhosts(VecFieldT& vec, int const levelNumber, double const fillTime)
    {
        auto schedule = findSchedule(vec.name(), levelNumber);
        if (schedule)
        {
            (*schedule)->fillData(fillTime);
        }
        else
        {
            throw std::runtime_error("no schedule for " + vec.name());
        }
    }




    void createInitSchedules(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
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
            refiner.add(algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy),
                        levelNumber);
        }
    }




    void initialize(int levelNumber, double initDataTime) const
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
};

} // namespace PHARE


#endif
