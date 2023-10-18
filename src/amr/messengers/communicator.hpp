#ifndef PHARE_QUANTITY_REFINER_HPP
#define PHARE_QUANTITY_REFINER_HPP


#include "amr/messengers/hybrid_messenger_info.hpp"

#include "amr/data/field/field_variable_fill_pattern.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"
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
    static constexpr auto is_refiner = std::is_same_v<ComType, Refiner>;

    template<typename ComType>
    static constexpr auto is_synchronizer = std::is_same_v<ComType, Synchronizer>;




    template<typename ComType, std::size_t dimension = 0>
    class Communicator
    {
    private:
        using Schedule  = typename ComType::schedule_type;
        using Algorithm = typename ComType::algorithm_type;
        std::map<int, std::map<Algorithm* const, std::shared_ptr<Schedule>>> schedules_;


    public:
        std::vector<std::unique_ptr<Algorithm>> algos;

        auto& add_algorithm()
        {
            if constexpr (is_refiner<ComType>)
                return algos.emplace_back(std::make_unique<SAMRAI::xfer::RefineAlgorithm>());

            if constexpr (is_synchronizer<ComType>)
                return algos.emplace_back(std::make_unique<SAMRAI::xfer::CoarsenAlgorithm>(
                    SAMRAI::tbox::Dimension{dimension}));
        }


        /**
         * @brief findSchedule returns the schedule at a given level number if there is one
         * (optional).
         */
        auto& findSchedule(std::unique_ptr<Algorithm> const& algo, int levelNumber) const
        {
            if (!schedules_.count(levelNumber))
                throw std::runtime_error("no schedule for level " + std::to_string(levelNumber));

            if (!schedules_.at(levelNumber).count(algo.get()))
                throw std::runtime_error("Algorithm has not been registered with Communicator");

            return schedules_.at(levelNumber).at(algo.get());
        }



        /**
         * @brief add is used to add a refine schedule for the given level number.
         * Note that already existing schedules at this level number are overwritten.
         */
        void add(std::unique_ptr<Algorithm> const& algo, std::shared_ptr<Schedule>& schedule,
                 int levelNumber)
        {
            // for shared border node value sync
            schedule->setDeterministicUnpackOrderingFlag(true);

            schedules_[levelNumber][algo.get()] = schedule;
        }

        void add(std::unique_ptr<Algorithm> const& algo, std::shared_ptr<Schedule>&& schedule,
                 int levelNumber)
        {
            add(algo, schedule, levelNumber);
        }
    };




    /**
     * @Brief This overload creates a Refiner for communication with both spatial and
     * time interpolation. Data is communicated from the model vector field defined at
     * time t=n+1 and its version at time t=n (oldModel), onto the `ghost` vector field.
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
     * @return the function returns a Refiner
     */
    template<typename ResourcesManager>
    Communicator<Refiner>
    makeRefiner(VecFieldDescriptor const& ghost, VecFieldDescriptor const& model,
                VecFieldDescriptor const& oldModel, std::shared_ptr<ResourcesManager> const& rm,
                std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp)
    {
        constexpr auto dimension = ResourcesManager::dimension;
        auto variableFillPattern = FieldFillPattern<dimension>::make_shared(refineOp);

        Communicator<Refiner> com;

        auto registerRefine
            = [&rm, &com, &refineOp, &timeOp](std::string const& ghost_, std::string const& model_,
                                              std::string const& oldModel_, auto& fillPattern) {
                  auto src_id  = rm->getID(ghost_);
                  auto dest_id = rm->getID(ghost_);
                  auto new_id  = rm->getID(model_);
                  auto old_id  = rm->getID(oldModel_);

                  if (src_id && dest_id && old_id)
                  {
                      com.add_algorithm()->registerRefine(
                          *dest_id, // dest
                          *src_id,  // source at same time
                          *old_id,  // source at past time (for time interp)
                          *new_id,  // source at future time (for time interp)
                          *dest_id, // scratch
                          refineOp, timeOp, fillPattern);
                  }
              };

        registerRefine(ghost.xName, model.xName, oldModel.xName, variableFillPattern);
        registerRefine(ghost.yName, model.yName, oldModel.yName, variableFillPattern);
        registerRefine(ghost.zName, model.zName, oldModel.zName, variableFillPattern);

        return com;
    }


    /**
     * @brief creates a Refiner for a scalar quantity with time refinement
     */
    template<typename ResourcesManager>
    Communicator<Refiner> makeRefiner(FieldDescriptor const& ghost, FieldDescriptor const& model,
                                      FieldDescriptor const& oldModel,
                                      std::shared_ptr<ResourcesManager> const& rm,
                                      std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                                      std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp)
    {
        constexpr auto dimension = ResourcesManager::dimension;
        auto variableFillPattern = FieldFillPattern<dimension>::make_shared(refineOp);

        Communicator<Refiner> com;

        auto registerRefine
            = [&rm, &com, &refineOp, &timeOp](std::string const& ghost_, std::string const& model_,
                                              std::string const& oldModel_, auto& fillPattern) {
                  auto src_id  = rm->getID(ghost_);
                  auto dest_id = rm->getID(ghost_);
                  auto new_id  = rm->getID(model_);
                  auto old_id  = rm->getID(oldModel_);

                  if (src_id && dest_id && old_id)
                  {
                      com.add_algorithm()->registerRefine(
                          *dest_id, // dest
                          *src_id,  // source at same time
                          *old_id,  // source at past time (for time interp)
                          *new_id,  // source at future time (for time interp)
                          *dest_id, // scratch
                          refineOp, timeOp, fillPattern);
                  }
              };

        registerRefine(ghost, model, oldModel, variableFillPattern);

        return com;
    }


    /**
     * @brief this overload creates a Refiner for communication without time interpolation
     * and from one quantity to the same quantity. It is typically used for initialization.
     */
    template<typename ResourcesManager>
    Communicator<Refiner> makeRefiner(VecFieldDescriptor const& descriptor,
                                      std::shared_ptr<ResourcesManager> const& rm,
                                      std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
    {
        return makeRefiner(descriptor, descriptor, rm, refineOp);
    }

    /**
     * @brief this overload creates a Refiner for communication without time interpolation
     * and from one quantity to another quantity.
     */
    template<typename ResourcesManager>
    Communicator<Refiner>
    makeRefiner(VecFieldDescriptor const& source, VecFieldDescriptor const& destination,
                std::shared_ptr<ResourcesManager> const& rm,
                std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, bool ghosts = false)
    {
        constexpr auto dimension = ResourcesManager::dimension;
        auto variableFillPattern = FieldFillPattern<dimension>::make_shared(refineOp);
        Communicator<Refiner> com;

        auto registerRefine
            = [&com, &rm, &refineOp, ghosts](std::string src, std::string dst, auto& fillPattern) {
                  auto idSrc  = rm->getID(src);
                  auto idDest = rm->getID(dst);
                  if (idSrc and idDest)
                  {
                      if (ghosts)
                          com.add_algorithm()->registerRefine(*idDest, *idSrc, *idDest, refineOp,
                                                              fillPattern);
                      else
                          com.add_algorithm()->registerRefine(*idDest, *idSrc, *idDest, refineOp);
                  }
              };
        registerRefine(source.xName, destination.xName, variableFillPattern);
        registerRefine(source.yName, destination.yName, variableFillPattern);
        registerRefine(source.zName, destination.zName, variableFillPattern);

        return com;
    }

    template<typename ResourcesManager>
    Communicator<Refiner> makeRefiner(std::string const& dest, std::string const& src,
                                      std::shared_ptr<ResourcesManager> const& rm,
                                      std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, bool)
    {
        Communicator<Refiner> com;

        auto idSrc  = rm->getID(src);
        auto idDest = rm->getID(dest);
        if (idSrc and idDest)
        {
            com.add_algorithm()->registerRefine(*idDest, *idSrc, *idDest, refineOp);
        }
        return com;
    }


    /**
     * @brief This overload of makeRefiner creates a Refiner for communication from one
     * scalar quantity to itself without time interpolation.
     */
    template<typename ResourcesManager>
    Communicator<Refiner> makeRefiner(std::string const& name,
                                      std::shared_ptr<ResourcesManager> const& rm,
                                      std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp)
    {
        return makeRefiner(name, name, rm, refineOp);
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

        auto registerCoarsen = [&com, &rm, &coarsenOp](std::string name) {
            auto id = rm->getID(name);
            if (id)
            {
                com.add_algorithm()->registerCoarsen(*id, *id, coarsenOp);
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
        {
            com.add_algorithm()->registerCoarsen(*id, *id, coarsenOp);
        }
        return com;
    }

} // namespace amr
} // namespace PHARE


#endif
