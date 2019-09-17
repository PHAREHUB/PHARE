#ifndef PHARE_COMMUNICATORS_H
#define PHARE_COMMUNICATORS_H

#include "quantity_communicator.h"

#include <SAMRAI/hier/RefineOperator.h>

#include <map>
#include <memory>
#include <string>



namespace PHARE
{
namespace amr_interface
{
    enum class RefinerType {
        GhostField,
        InitField,
        InitInteriorPart,
        LevelBorderParticles,
        InteriorGhostParticles
    };


    /**
     * @brief The Communicators class is used by a Messenger to manipulate SAMRAI algorithms and
     * schedules It contains a QuantityCommunicator for all quantities registered to the Messenger
     * for ghost, init etc.
     */
    template<RefinerType Type>
    class RefinerPool
    {
    public:
        /**
         * @brief add a QuantityCommunicator to the communicators based on the given Descriptor and
         * the refinement operator refineOp. The method uses the ResourcesManager to check the
         * quantities in the Descriptor are properly registered. The created QuantityCommunicator is
         * associated with the given key.
         *
         * This overload is used for communications that do not use time interpolation
         */
        template<typename Descriptor, typename ResourcesManager>
        void add(Descriptor const& descriptor,
                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, std::string key,
                 std::shared_ptr<ResourcesManager> const& rm)
        {
            auto const [it, success]
                = refiners_.insert({key, makeRefiner(descriptor, rm, refineOp)});

            if (!success)
                throw std::runtime_error(key + " is already registered");
        }


        /**
         * @brief add a QuantityCommunicator to the Communicators based on the given descroptors and
         * spatial refinement operator and time interpolation operator. The created
         * QuantityCommunicator is associated with the given key.
         *
         * This overload is used for communications that involve time interpolation. In PHARE those
         * are only for ghost nodes, which is only for VecField E and B.
         */
        template<typename ResourcesManager>
        void add(VecFieldDescriptor const& ghostDescriptor,
                 VecFieldDescriptor const& modelDescriptor,
                 VecFieldDescriptor const& oldModelDescriptor,
                 std::shared_ptr<ResourcesManager> const& rm,
                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                 std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> timeOp, std::string key)
        {
            auto const [it, success]
                = refiners_.insert({key, makeRefiner(ghostDescriptor, modelDescriptor,
                                                     oldModelDescriptor, rm, refineOp, timeOp)});
            if (!success)
                throw std::runtime_error(key + " is already registered");
        }




        /**
         * @brief registerLevel registers a level of the hierarchy to all QuantityCommunicators in
         * the Communicators.
         *
         * For each QuantityCommunicator, the method takes the RefineAlgorithm and creates a
         * schedule by calling one of the createSchedule() overloads. The specific overload that is
         * called depends on the (compile-time) nature of the Communicators.
         *
         */
        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
        {
            for (auto& [_, refiner] : refiners_)
            {
                auto& algo       = refiner.algo;
                auto levelNumber = level->getLevelNumber();


                // for GhostField we need schedules that take on the level where there is an overlap
                // (there is always for patches lying inside the level)
                // and goes to coarser level where there is not (patch lying on the level border)
                if constexpr (Type == RefinerType::GhostField)
                {
                    auto schedule = algo->createSchedule(
                        level, level->getNextCoarserHierarchyLevelNumber(), hierarchy);
                    refiner.add(schedule, levelNumber);
                }

                // this createSchedule overload is used to initialize fields.
                // note that here we must take that createsSchedule() overload and put nullptr as
                // src since we want to take from coarser level everywhere. using the createSchedule
                // overload that takes level, next_coarser_level only would result in interior ghost
                // nodes to be filled with interior of neighbor patches but there is nothing there.
                else if constexpr (Type == RefinerType::InitField)
                {
                    refiner.add(algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy),
                                levelNumber);
                }


                // here we create the schedule that will intialize the particles that lie within the
                // interior of the patches (no ghost, no coarse to fine). We take almost the same
                // overload as for fields above but the version that takes a PatchLevelFillPattern.
                // Here the PatchLevelInteriorFillPattern is used because we want to fill particles
                // only within the interior of the patches of the level. The reason is that filling
                // the their ghost regions with refined particles would not ensure the ghosts to be
                // clones of neighbor patches particles if the splitting from coarser levels is not
                // deterministic.
                else if constexpr (Type == RefinerType::InitInteriorPart)
                {
                    refiner.add(algo->createSchedule(
                                    std::make_shared<SAMRAI::xfer::PatchLevelInteriorFillPattern>(),
                                    level, nullptr, levelNumber - 1, hierarchy),
                                levelNumber);
                }

                // here we create a schedule that will refine particles from coarser level and put
                // them into the level coarse to fine boundary. These are the levelGhostParticlesOld
                // particles. we thus take the same createSchedule overload as above but pass it a
                // PatchLevelBorderFillPattern.
                else if constexpr (Type == RefinerType::LevelBorderParticles)
                {
                    refiner.add(algo->createSchedule(
                                    std::make_shared<SAMRAI::xfer::PatchLevelBorderFillPattern>(),
                                    level, nullptr, levelNumber - 1, hierarchy),
                                levelNumber);
                }

                // this branch is used to create a schedule that will transfer particles into the
                // patches' ghost zones.
                else if constexpr (Type == RefinerType::InteriorGhostParticles)
                {
                    refiner.add(algo->createSchedule(level), levelNumber);
                }
            }
        }




        /**
         * @brief initialize is used to initialize data on the level for all quantities in the
         * Communicators.
         *
         * Basically the method just takes all QuantityCommunicator one by one, find the schedule
         * associated with the given level number and executes fillData().
         *
         * The method registerLevel must have been called before for the given levelNumber otherwise
         * no schedule will be found
         */
        void fill(int levelNumber, double initDataTime) const
        {
            for (auto& [key, communicator] : refiners_)
            {
                if (communicator.algo == nullptr)
                {
                    throw std::runtime_error("Algorithm is nullptr");
                }

                auto schedule = communicator.findSchedule(levelNumber);
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




        /**
         * @brief regrid is used to execute a regridding schedule for all quantities in the pool.
         */
        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            const int levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            double const initDataTime)
        {
            for (auto& [key, refiner] : refiners_)
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




        /**
         * @brief filldGhosts is used to fill the ghost nodes of all the given VecField.
         *
         * The VecField must have been registered before to the pool for ghost filling, and the
         * method createGhostSchedule must have been called for the given levelNumber otherwise no
         * refine schedule will be found and the method will throw an axception.
         */
        template<typename VecFieldT>
        void fill(VecFieldT& vec, int const levelNumber, double const fillTime)
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
        findSchedule_(std::string const& name, int levelNumber)
        {
            if (auto mapIter = refiners_.find(name); mapIter != std::end(refiners_))
            {
                return mapIter->second.findSchedule(levelNumber);
            }
            else
            {
                return std::nullopt;
            }
        }



        std::map<std::string, Communicator<Refiner>> refiners_;
    };




    template<std::size_t dimension>
    class SynchronizerPool
    {
    public:
        template<typename Descriptor, typename ResourcesManager>
        void add(Descriptor const& descriptor,
                 std::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsenOp, std::string key,
                 std::shared_ptr<ResourcesManager> const& rm)
        {
            auto const [it, success] = synchronizers_.insert(
                {key, makeSynchronizer<ResourcesManager, dimension>(descriptor, rm, coarsenOp)});
        }




        template<typename ResourcesManager>
        void add(VecFieldDescriptor const& descriptor, std::shared_ptr<ResourcesManager> const& rm,
                 std::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsenOp, std::string key)
        {
            auto const [it, success] = synchronizers_.insert(
                {key, makeSynchronizer<ResourcesManager, dimension>(descriptor, rm, coarsenOp)});

            if (!success)
                throw std::runtime_error(key + " is already registered");
        }



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
        {
            for (auto& [_, synchronizer] : synchronizers_)
            {
                auto& algo              = synchronizer.algo;
                auto levelNumber        = level->getLevelNumber();
                auto const& coarseLevel = hierarchy->getPatchLevel(levelNumber - 1);

                auto schedule = algo->createSchedule(coarseLevel, level);
                synchronizer.add(std::move(schedule), levelNumber);
            }
        }



        void sync(int const levelNumber)
        {
            for (auto& [key, synchronizer] : synchronizers_)
            {
                if (synchronizer.algo == nullptr)
                {
                    throw std::runtime_error("Algorithm is nullptr");
                }

                auto schedule = synchronizer.findSchedule(levelNumber);
                if (schedule)
                {
                    (*schedule)->coarsenData();
                }
                else
                {
                    throw std::runtime_error("Error - schedule cannot be found for this level");
                }
            }
        }

    private:
        std::map<std::string, Communicator<Synchronizer, dimension>> synchronizers_;
    };


} // namespace amr_interface
} // namespace PHARE

#endif
