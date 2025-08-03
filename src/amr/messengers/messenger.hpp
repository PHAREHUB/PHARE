
#ifndef PHARE_MESSENGER_HPP
#define PHARE_MESSENGER_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include "messenger_info.hpp"


namespace PHARE
{
namespace amr
{
    /**
     * @brief The IMessenger class is an interface used to communicate data between patches and
     * levels. It is used by the MultiPhysicsIntegrator and the ISolver objects. These
     * communications are of several categories:
     *
     * - initializing a level from a coarser one (done from the MultiPhysicsIntegrator)
     * - regridding, i.e. initializing a level from a former one at the same refinement level and a
     * coarser one (done from the MultiPhysicsIntegrator)
     * - filling ghost nodes from the coarser level (done by concrete ISolver objects)
     * - synchronizing quantities from the current level to the coarser one. (done by the
     * MultiPhysicsIntegrator)
     *
     * These kinds of communications are general to the evolution of a plasma hierarchy, no matter
     * the nature of the models at the fine and coarse levels in between which the messengers take
     * place.
     *
     * However, for each of these operations several quantities are concerned and their nature
     * depends on the nature of the models at the fine and coarse levels. For instance, messengers
     * of data between two MHD levels will not involve the same quantities for each of the above
     * cases as would, say, messengers of data beween two hybrid levels.
     *
     * Therefore, specializing the IMessenger interface with a concrete implementation is required
     * for each combination of (fine model, coarse model) involved in the PatchHierarchy.
     *
     * Concrete specializations are:
     * - HybridMessenger
     * - MHDMessenger
     *
     * IMessenger and subclasses exist also to hide all the SAMRAI details involved in communicating
     * data in a PatchHierarchy and to provide a simple interface to the MultiPhysicsIntegrator and
     * ISolver.
     *
     */
    template<typename IPhysicalModel>
    class IMessenger
    {
    public:
        virtual std::string name() = 0;

        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() = 0;
        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner()   = 0;


        /**
         * @brief allocate is used to allocate resources used internally by the concrete IMessenger.
         *
         * This routine will be called by the MultiPhysicsIntegrator while in
         * initializeLevelData. Internal resources owned by an IMessengers can be for instance
         * temporal temporary versions of model variable used for time interpolations. Other kind of
         * internal resources can be temporary variables used in multi-model messengers involving
         * quantities of different nature. Any internal data used by the messenger that needs to be
         * allocated onto Patches needs to be allocated in this methods.
         *
         * @param patch is the samrai Patch onto which we want to allocate the data
         * @param allocateTime is the time at which data is allocated
         */
        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const = 0;


        /**
         * @brief registerQuantities is used to tell the messenger which data will be communicated
         * between coarser and finer levels. This method is called only once when the Messenger
         * is created. Only information regarding which data will be communicated is given, no
         * information about which levels are concerned because that one needs to be updated
         * regularly as the patch hierarchy evolves.
         *
         * registerQuantities needs information from both the coarse and fine levels involved in the
         * Messenger. The method takes ownership of abstract IMessengerInfo structures. Concrete
         * implementation of setup, e.g. HybridMessenger will need to cast these into their derived
         * class of IMessengerInfo, e.g. HybridMessengerInfo.
         *
         * @param fromCoarserInfo is an abstract object of which concrete derived structure
         * encapsulates the names of all data involved in this messenger that comes from the coarse
         * level.
         *
         * @param fromFinerInfo is an abstract object of which concrete derived structure
         * encapsulates the names of all data involved in this messenger that comes from the fine
         * level.
         *
         * precondition: the variable names contained in concrete IMessengerInfo structures must
         * correspond to the names of model and solver variables that have already been registered.
         * In case one name does not correspond to a known registered variable, the code will throw
         * an exception.
         */
        virtual void registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                        std::unique_ptr<IMessengerInfo> fromFinerInfo)
            = 0;



        /**
         * @brief registerLevel tells the messenger which level in the hierarchy data will be
         * communicated to. When the messenger involves communication from coarser grid, it is from
         * levelNumber-1.
         *
         * This method needs to be called whenever initializing a new level in the
         * hierarchy so that a valid Messenger be used at this level. The method registerQuantities
         * must be called before this method. The levelnumber needs to be the one of a level in the
         * given hierarchy.
         */
        virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   int const levelNumber)
            = 0;



        /**
         * @brief regrid performs communications necessary to fill a new level when regridding.
         * After this method is called, all variables needed on this new level are filled where they
         * need to be. Whether they belong to the model, the solver or the messenger itself.
         *
         * @param hierarchy is the patch hierarchy being regridded
         * @param levelNumber is the number of the new level
         * @param oldLevel is the old level being removed from the hierarchy and from which we can
         * take data where the new level overlaps it
         * @param initDataTime is the time of the regridding
         */
        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            int const levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            IPhysicalModel& model, double const initDataTime)
            = 0;




        /**
         * @brief initLevel is used by the MultiPhysicsIntegrator to initialize a level at time
         * initDataTime from a coarser level iLevel-1. The initialization is performed on all
         * quantities registered during registerQuantities().
         */
        virtual void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               double const initDataTime)
            = 0;




        /**
         * @brief this method is used by the MultiPhysicsIntegrator at the first step of a
         * substepping cycle because some concrete IMessenger object may need to perform some
         * operations at this time.
         *
         * Typically a concrete IMessenger may need to get coarser level data only at the first step
         * and then use it for time interpolation at each substep.
         *
         * @param model
         * @param level
         * @param time
         */
        virtual void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               double const currentTime, double const prevCoarserTime,
                               double const newCoarserTime)
            = 0;




        /**
         * @brief lastStep is a method called by the MultiPhysicsIntegrator at the last step of the
         * substepping cycle. Concrete IMessenger objects may need to perform operations at this
         * step.
         *
         * @param model
         */
        virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) = 0;




        /**
         * @brief prepareStep is used by the MultiPhysicsIntegrator just before calling the solver
         * advanceLevel function. The method sets the Messenger in a state that is ready for the
         * solver to actually use it during advanceLevel.
         */
        virtual void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                 double currentTime)
            = 0;



        virtual void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                    double const initDataTime)
            = 0;



        virtual void synchronize(SAMRAI::hier::PatchLevel& level) = 0;

        virtual void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                            double const syncTime)
            = 0;

        virtual void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                     double const time)
            = 0;


        /**
         * @brief fineModelName returns the name of the fine model involved in the messenger
         * @return
         */
        virtual std::string fineModelName() const = 0;



        /**
         * @brief coarseModelName returns the name of the coarse model involved i the messenger
         * @return
         */
        virtual std::string coarseModelName() const = 0;

        virtual ~IMessenger() = default;
    };


    template<typename IPhysicalModel>
    bool areCompatible(IMessenger<IPhysicalModel> const& messenger, IPhysicalModel const& model)
    {
        return model.name() == messenger.fineModelName()
               || model.name() == messenger.coarseModelName();
    }

} // namespace amr
} // namespace PHARE
#endif
