#ifndef PHARE_QUANTITY_REFINER_H
#define PHARE_QUANTITY_REFINER_H


#include "evolution/messengers/hybrid_messenger_info.h"

#include "data/field/refine/field_refine_operator.h"
#include "data/field/time_interpolate/field_linear_time_interpolate.h"
#include "data/particles/refine/particles_data_split.h"
#include <SAMRAI/xfer/PatchLevelBorderFillPattern.h>
#include <SAMRAI/xfer/PatchLevelInteriorFillPattern.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>

#include <map>
#include <memory>
#include <optional>

namespace PHARE
{
namespace amr_interface
{
    /**
     * @brief The Refiner struct encapsulate a SAMRAI algorithm for a quantity and all the schedules
     * associated to the levels in the hierarchy.
     */
    class QuantityCommunicator
    {
    public:
        QuantityCommunicator()
            : algo{std::make_unique<SAMRAI::xfer::RefineAlgorithm>()}
        {
        }


        /**
         * @brief findSchedule returns the schedule at a given level number if there is one
         * (optional).
         */
        std::optional<std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
        findSchedule(int levelNumber) const
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
     * passing it the IDs of the ghost, model and old model patch datas associated to each component
     * of the vector field.
     *
     *
     * @param ghost is the VecFieldDescriptor of the VecField that needs its ghost nodes filled
     * @param model is the VecFieldDescriptor of the model VecField from which data is taken (at
     * time t_coarse+dt_coarse)
     * @param oldModel is the VecFieldDescriptor of the model VecField from which data is taken at
     * time t_coarse
     * @param rm is the ResourcesManager
     * @param refineOp is the spatial refinement operator
     * @param timeOp is the time interpolator
     *
     * @return the function returns a QuantityRefiner which may be stored in a RefinerPool and to
     * which later schedules will be added.
     */
    template<typename ResourcesManager>
    QuantityCommunicator
    makeCommunicator(VecFieldDescriptor const& ghost, VecFieldDescriptor const& model,
                     VecFieldDescriptor const& oldModel, ResourcesManager const& rm,
                     std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                     std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp)
    {
        QuantityCommunicator refiner;

        auto registerRefine = [&rm, &refiner, &refineOp, &timeOp](std::string const& modelName,
                                                                  std::string const& ghostName,
                                                                  std::string const& oldModelName) {
            auto src_id  = rm->getID(modelName);
            auto dest_id = rm->getID(ghostName);
            auto old_id  = rm->getID(oldModelName);

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




    /**
     * @brief makeInitRefiner is similar to makeGhostRefiner except the registerRefine() that is
     * called is the one that allows initialization of a vector field quantity.
     */
    template<typename ResourcesManager>
    QuantityCommunicator makeCommunicator(VecFieldDescriptor const& _name,
                                          ResourcesManager const& rm,
                                          std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
    {
        QuantityCommunicator refiner;

        auto registerRefine = [&refiner, &rm, &refineOp](std::string name) //
        {
            auto id = rm->getID(name);
            if (id)
            {
                refiner.algo->registerRefine(*id, *id, *id, refineOp);
            }
        };

        registerRefine(_name.xName);
        registerRefine(_name.yName);
        registerRefine(_name.zName);

        return refiner;
    }




    /**
     * @brief makeInitRefiner is similar to makeGhostRefiner except the registerRefine() that is
     * called is the one that allows initialization of a field quantity.
     */
    template<typename ResourcesManager>
    QuantityCommunicator makeCommunicator(std::string const& name, ResourcesManager const& rm,
                                          std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
    {
        QuantityCommunicator communicator;

        auto id = rm->getID(name);
        if (id)
        {
            communicator.algo->registerRefine(*id, *id, *id, refineOp);
        }

        return communicator;
    }



} // namespace amr_interface
} // namespace PHARE


#endif
