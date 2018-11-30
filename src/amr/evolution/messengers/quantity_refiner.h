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
 * @brief The Refiner struct encapsulate a SAMRAI algorithm for a quantity and all the schedules
 * associated to the levels in the hierarchy.
 */
class QuantityRefiner
{
public:
    QuantityRefiner()
        : algo{std::make_unique<SAMRAI::xfer::RefineAlgorithm>()}
    {
    }


    /**
     * @brief findSchedule returns the schedule at a given level number if there is one (optional).
     */
    std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>> findSchedule(int levelNumber) const
    {
        if (auto mapIter = schedules_.find(levelNumber); mapIter != std::end(schedules_))
        {
            return mapIter->second;
        }
        else
        {
            return std::nullopt;
        }
    }



    /**
     * @brief add is used to add a refine schedule for the given level number.
     * Note that already existing schedules at this level number are overwritten.
     */
    void add(std::shared_ptr<SAMRAI::xfer::RefineSchedule> schedule, int levelNumber)
    {
        schedules_[levelNumber] = std::move(schedule);
    }



    std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> algo;

private:
    std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> schedules_;
};




/**
 * @brief makeGhostRefiner creates a QuantityRefiner for ghost filling of a VecField.
 *
 * The method basically calls registerRefine() on the QuantityRefiner algorithm,
 * passing it the IDs of the ghost, model and old model patch datas associated to each component of
 * the vector field.
 *
 *
 * @param ghost is the VecFieldDescriptor of the VecField that needs its ghost nodes filled
 * @param model is the VecFieldDescriptor of the model VecField from which data is taken (at time
 * t_coarse+dt_coarse)
 * @param oldModel is the VecFieldDescriptor of the model VecField from which data is taken at time
 * t_coarse
 * @param rm is the ResourcesManager
 * @param refineOp is the spatial refinement operator
 * @param timeOp is the time interpolator
 *
 * @return the function returns a QuantityRefiner which may be stored in a RefinerPool and to which
 * later schedules will be added.
 */
template<typename ResourcesManager>
QuantityRefiner makeGhostRefiner(VecFieldDescriptor const& ghost, VecFieldDescriptor const& model,
                                 VecFieldDescriptor const& oldModel, ResourcesManager const& rm,
                                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                                 std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp)
{
    QuantityRefiner refiner;

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
/**
 * @brief makeInitRefiner is similar to makeGhostRefiner except the registerRefine() that is called
 * is the one that allows initialization of a vector field quantity.
 */
QuantityRefiner makeInitRefiner(VecFieldDescriptor const& name, ResourcesManager const& rm,
                                std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
{
    QuantityRefiner refiner;

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




/**
 * @brief The RefinerPool class is used by a Messenger to manipulate SAMRAI algorithms and schedules
 * It contains a QuantityRefiner for all quantities registered to the Messenger for ghost, init etc.
 */
class RefinerPool
{
public:
    /**
     * @brief add is used to add a QuantityRefiner to the pool, associated with the key qtyName.
     */
    void add(QuantityRefiner&& qtyRefiner, std::string const& qtyName);



    /**
     * @brief createGhostSchedulesis adds a ghost filling schedule to all QuantityRefiner of the
     * pool. Once this method is called, all quantities in the pool can have their ghost filled on
     * the given level using fillVecFieldGhosts().
     */
    void createGhostSchedules(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              std::shared_ptr<SAMRAI::hier::PatchLevel>& level);



    /**
     * @brief createInitSchedules adds initialization schedule to all quantities in the pool on the
     * given level. This method needs to be called to later use initialize() on the given level.
     */
    void createInitSchedules(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             std::shared_ptr<SAMRAI::hier::PatchLevel> const& level);




    /**
     * @brief initialize is used to initialize all quantities registered in the pool.
     *
     * Basically the method just takes all QuantityRefiner one by one, find the schedule associated
     * with the given level number and executes fillData().
     *
     * The method createInitSchedule() must have been called before for that level otherwise
     * initialize() will not find schedules and will throw a run time exception.
     */
    void initialize(int levelNumber, double initDataTime) const;




    /**
     * @brief regrid is used to execute a regridding schedule for all quantities in the pool.
     */
    virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        const int levelNumber,
                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                        double const initDataTime);




    /**
     * @brief fillVecFieldGhosts is used to fill the ghost nodes of all the given VecField.
     *
     * The VecField must have been registered before to the pool for ghost filling, and the method
     * createGhostSchedule must have been called for the given levelNumber otherwise no refine
     * schedule will be found and the method will throw an axception.
     */
    template<typename VecFieldT>
    void fillVecFieldGhosts(VecFieldT& vec, int const levelNumber, double const fillTime)
    {
        auto schedule = findSchedule_(vec.name(), levelNumber);
        if (schedule)
        {
            (*schedule)->fillData(fillTime);
        }
        else
        {
            throw std::runtime_error("no schedule for " + vec.name());
        }
    }



private:
    std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
    findSchedule_(std::string const& name, int levelNumber);



    std::map<std::string, QuantityRefiner> qtyRefiners_;
};

} // namespace PHARE


#endif
