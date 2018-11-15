
#ifndef PHARE_TRANSACTION_H
#define PHARE_TRANSACTION_H

#include <string>

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include "physical_models/physical_model.h"
#include "transaction_info.h"


namespace PHARE
{
/**
 * @brief The ITransaction class is an interface used to communicate data between patches and
 * levels. It is used by the MultiPhysicsIntegrator and the ISolver objects. These communications
 * are of several categories:
 *
 * - initializing a level from a coarser one (done from the MultiPhysicsIntegrator)
 * - regridding, i.e. initializing a level from a former one at the same refinement level and a
 * coarser one (done from the MultiPhysicsIntegrator)
 * - filling ghost nodes from the coarser level (done by concrete ISolver objects)
 * - synchronizing quantities from the current level to the coarser one. (done by the
 * MultiPhysicsIntegrator)
 *
 * These kinds of communications are general to the evolution of a plasma hierarchy, no matter the
 * nature of the models at the fine and coarse levels in between which the transactions take place.
 *
 * However, for each of these operations several quantities are concerned and their nature depends
 * on the nature of the models at the fine and coarse levels. For instance, transactions of data
 * between two MHD levels will not involve the same quantities for each of the above
 * cases as would, say, transactions of data beween two hybrid levels.
 *
 * Therefore, specializing the ITransaction interface with a concrete implementation is required for
 * each combination of (fine model, coarse model) involved in the PatchHierarchy.
 *
 * Concrete specializations are:
 * - HybridTransaction
 * - MHDTransaction
 *
 * ITransaction and subclasses exist also to hide all the SAMRAI details involved in communicating
 * data in a PatchHierarchy and to provide a simple interface to the MultiPhysicsIntegrator and
 * ISolver.
 *
 */
class ITransaction
{
public:
    virtual std::string name() = 0;

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() = 0;
    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner()   = 0;


    /**
     * @brief allocate is used to allocate resources used internally by the concrete ITransaction.
     *
     * This routine will be called by the MultiPhysicsIntegrator while in
     * initializeLevelData. Internal resources owned by an ITransactions can be for instance
     * temporal temporary versions of model variable used for time interpolations. Other kind of
     * internal resources can be temporary variables used in multi-model transactions involving
     * quantities of different nature. Any internal data used by the transaction that needs to be
     * allocated onto Patches needs to be allocated in this methods.
     *
     * @param patch is the samrai Patch onto which we want to allocate the data
     * @param allocateTime is the time at which data is allocated
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const = 0;


    /**
     * @brief registerQuantities is used to tell the transaction which data will be communicated
     * between coarser and finer levels. This method is called only once when the Transaction
     * is created. Only information regarding which data will be communicated is given, no
     * information about which levels are concerned because that one needs to be updated
     * regularly as the patch hierarchy evolves.
     *
     * registerQuantities needs information from both the coarse and fine levels involved in the
     * Transaction. The method takes ownership of abstract ITransactionInfo structures. Concrete
     * implementation of setup, e.g. HybridTransaction will need to cast these into their derived
     * class of ITransactionInfo, e.g. HybridTransactionInfo.
     *
     * @param fromCoarserInfo is an abstract object of which concrete derived structure encapsulates
     * the names of all data involved in this transaction that comes from the coarse level.
     *
     * @param fromFinerInfo is an abstract object of which concrete derived structure encapsulates
     * the names of all data involved in this transaction that comes from the fine level.
     *
     * precondition: the variable names contained in concrete ITransactionInfo structures must
     * correspond to the names of model and solver variables that have already been registered.
     * In case one name does not correspond to a known registered variable, the code will throw an
     * exception.
     */
    virtual void registerQuantities(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                                    std::unique_ptr<ITransactionInfo> fromFinerInfo)
        = 0;



    /**
     * @brief registerLevel tells the transaction which level in the hierarchy data will be
     * communicated to. When the transaction involves communication from coarser grid, it is from
     * levelNumber-1.
     *
     * This method needs to be called whenever initializing a new level in the
     * hierarchy so that a valid Transaction be used at this level. The method registerQuantities
     * must be called before this method. The levelnumber needs to be the one of a level in the
     * given hierarchy.
     */
    virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               int const levelNumber)
        = 0;



    /**
     * @brief regrid performs communications necessary to fill a new level when regridding.
     * After this method is called, all variables needed on this new level are filled where they
     * need to be. Whether they belong to the model, the solver or the transaction itself.
     *
     * @param hierarchy is the patch hierarchy being regridded
     * @param levelNumber is the number of the new level
     * @param oldLevel is the old level being removed from the hierarchy and from which we can take
     * data where the new level overlaps it
     * @param initDataTime is the time of the regridding
     */
    virtual void
    regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy, const int levelNumber,
           std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel, double const initDataTime)
        = 0;




    /**
     * @brief initLevel is used by the MultiPhysicsIntegrator to initialize a level at time
     * initDataTime from a coarser level iLevel-1. The initialization is performed on all quantities
     * registered during registerQuantities().
     *
     * @param iLevel
     * @param initDataTime
     */
    virtual void initLevel(int const iLevel, double const initDataTime) const = 0;




    /**
     * @brief this method is used by the MultiPhysicsIntegrator at the first step of a substepping
     * cycle because some concrete ITransaction object may need to perform some operations at this
     * time.
     *
     * Typically a concrete ITransaction may need to get coarser level data only at the first step
     * and then use it for time interpolation at each substep.
     *
     * @param model
     */
    virtual void firstStep(IPhysicalModel const& model) = 0;




    /**
     * @brief lastStep is a method called by the MultiPhysicsIntegrator at the last step of the
     * substepping cycle. Concrete ITransaction objects may need to perform operations at this step.
     *
     * @param model
     */
    virtual void lastStep(IPhysicalModel const& model) = 0;



    /**
     * @brief fineModelName returns the name of the fine model involved in the transaction
     * @return
     */
    virtual std::string fineModelName() const = 0;



    /**
     * @brief coarseModelName returns the name of the coarse model involved i the transaction
     * @return
     */
    virtual std::string coarseModelName() const = 0;



    virtual ~ITransaction() = default;
};



bool areCompatible(ITransaction const& transaction, IPhysicalModel const& model);



} // namespace PHARE
#endif
