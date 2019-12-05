#ifndef PHARE_QUANTITY_REFINER_H
#define PHARE_QUANTITY_REFINER_H


#include "amr/messengers/hybrid_messenger_info.h"

#include "amr/data/field/coarsening/field_coarsen_operator.h"
#include "amr/data/field/refine/field_refine_operator.h"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.h"
#include "amr/data/particles/refine/particles_data_split.h"
#include <SAMRAI/tbox/Dimension.h>
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#include <SAMRAI/xfer/CoarsenSchedule.h>
#include <SAMRAI/xfer/PatchLevelBorderFillPattern.h>
#include <SAMRAI/xfer/PatchLevelInteriorFillPattern.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"

#include <map>
#include <memory>
#include <optional>



namespace PHARE
{
namespace amr
{
    /* When determining which DataFactory to use to create the correct geometry
     *  Samrai compares existing items within the RefinerPool
     *   See: SAMRAI::xfer::RefineClasses::getEquivalenceClassIndex"
     *   If we do not specify a separate VariableFillPattern for each item sent to
     *    Algo->registerRefine : we run the risk of all items only using the first registered
     *    DataFactory, thus all items will receive the Geometry of that item.
     *    We "hack" this as there is a typeid()== check on the variablefill pattern types
     */
    class XVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
    {
    };

    class YVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
    {
    };

    class ZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
    {
    };

    struct Refiner
    {
        using schedule_type  = SAMRAI::xfer::RefineSchedule;
        using algorithm_type = SAMRAI::xfer::RefineAlgorithm;
    };

    struct Synchronizer
    {
        using schedule_type  = SAMRAI::xfer::CoarsenSchedule;
        using algorithm_type = SAMRAI::xfer::CoarsenAlgorithm;
    };

    template<typename ComType>
    using is_refiner = std::enable_if_t<std::is_same_v<ComType, Refiner>>;

    template<typename ComType>
    using is_synchronizer = std::enable_if_t<std::is_same_v<ComType, Synchronizer>>;




    template<typename ComType, std::size_t dimension = 0>
    class Communicator
    {
    private:
        using Schedule  = typename ComType::schedule_type;
        using Algorithm = typename ComType::algorithm_type;
        std::map<int, std::shared_ptr<Schedule>> schedules_;

    public:
        std::unique_ptr<Algorithm> algo;



        template<typename U = ComType, typename = is_refiner<U>>
        Communicator()
            : algo{std::make_unique<SAMRAI::xfer::RefineAlgorithm>()}
        {
        }



        template<typename Dummy = ComType, typename = is_synchronizer<Dummy>, typename = Dummy>
        Communicator()
            : algo{std::make_unique<SAMRAI::xfer::CoarsenAlgorithm>(
                SAMRAI::tbox::Dimension{dimension})}

        {
        }



        /**
         * @brief findSchedule returns the schedule at a given level number if there is one
         * (optional).
         */
        std::optional<std::shared_ptr<Schedule>> findSchedule(int levelNumber) const
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
        void add(std::shared_ptr<Schedule> schedule, int levelNumber)
        {
            schedules_[levelNumber] = std::move(schedule);
        }
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
    Communicator<Refiner>
    makeRefiner(VecFieldDescriptor const& ghost, VecFieldDescriptor const& model,
                VecFieldDescriptor const& oldModel, std::shared_ptr<ResourcesManager> const& rm,
                std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp)
    {
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> xVariableFillPattern
            = std::make_shared<XVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> yVariableFillPattern
            = std::make_shared<YVariableFillPattern>();
        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> zVariableFillPattern
            = std::make_shared<ZVariableFillPattern>();

        Communicator<Refiner> com;

        auto registerRefine
            = [&rm, &com, &refineOp, &timeOp](std::string const& model_, std::string const& ghost_,
                                              std::string const& oldModel_, auto& fillPattern) {
                  auto src_id  = rm->getID(model_);
                  auto dest_id = rm->getID(ghost_);
                  auto old_id  = rm->getID(oldModel_);

                  if (src_id && dest_id && old_id)
                  {
                      // dest, src, old, new, scratch
                      com.algo->registerRefine(*dest_id, // dest
                                               *src_id,  // source at same time
                                               *old_id,  // source at past time (for time interp)
                                               *src_id,  // source at future time (for time interp)
                                               *dest_id, // scratch
                                               refineOp, timeOp, fillPattern);
                  }
              };

        // register refine operators for each component of the vecfield
        registerRefine(ghost.xName, model.xName, oldModel.xName, xVariableFillPattern);
        registerRefine(ghost.yName, model.yName, oldModel.yName, yVariableFillPattern);
        registerRefine(ghost.zName, model.zName, oldModel.zName, zVariableFillPattern);

        return com;
    }




    /**
     * @brief makeInitRefiner is similar to makeGhostRefiner except the registerRefine() that is
     * called is the one that allows initialization of a vector field quantity.
     */
    template<typename ResourcesManager>
    Communicator<Refiner> makeRefiner(VecFieldDescriptor const& descriptor,
                                      std::shared_ptr<ResourcesManager> const& rm,
                                      std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
    {
        Communicator<Refiner> com;

        auto registerRefine = [&com, &rm, &refineOp](std::string name) //
        {
            auto id = rm->getID(name);
            if (id)
            {
                com.algo->registerRefine(*id, *id, *id, refineOp);
            }
        };
        registerRefine(descriptor.xName);
        registerRefine(descriptor.yName);
        registerRefine(descriptor.zName);

        return com;
    }




    /**
     * @brief makeInitRefiner is similar to makeGhostRefiner except the registerRefine() that is
     * called is the one that allows initialization of a field quantity.
     */
    template<typename ResourcesManager>
    Communicator<Refiner> makeRefiner(std::string const& name,
                                      std::shared_ptr<ResourcesManager> const& rm,
                                      std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
    {
        Communicator<Refiner> communicator;

        auto id = rm->getID(name);
        if (id)
        {
            communicator.algo->registerRefine(*id, *id, *id, refineOp);
        }

        return communicator;
    }



    /**
     * @brief makeInitRefiner is similar to makeGhostRefiner except the registerRefine() that is
     * called is the one that allows initialization of a vector field quantity.
     */
    template<typename ResourcesManager, std::size_t dimension>
    Communicator<Synchronizer, dimension>
    makeSynchronizer(VecFieldDescriptor const& descriptor,
                     std::shared_ptr<ResourcesManager> const& rm,
                     std::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsenOp)
    {
        Communicator<Synchronizer, dimension> com;

        auto registerCoarsen = [&com, &rm, &coarsenOp](std::string name) //
        {
            auto id = rm->getID(name);
            if (id)
            {
                com.algo->registerCoarsen(*id, *id, coarsenOp);
            }
        };

        registerCoarsen(descriptor.xName);
        registerCoarsen(descriptor.yName);
        registerCoarsen(descriptor.zName);

        return com;
    }




    template<typename ResourcesManager, std::size_t dimension>
    Communicator<Synchronizer, dimension>
    makeSynchronizer(std::string const& name, std::shared_ptr<ResourcesManager> const& rm,
                     std::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsenOp)
    {
        Communicator<Synchronizer, dimension> com;

        auto id = rm->getID(name);
        if (id)
            com.algo->registerCoarsen(*id, *id, coarsenOp);

        return com;
    }

} // namespace amr
} // namespace PHARE


#endif
