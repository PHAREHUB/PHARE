
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
 * @brief The ITransaction class is an interface used to communicate data between patches and level.
 * It is used by the MultiPhysicsIntegrator and the ISolver objects. ITransaction ...
 *
 * By hiding all the details related to SAMRAI algorithms and schedules, and by providing a simple
 * interface, ITransaction makes possible to perform communications in the MultiPhysicsIntegrator or
 * ISolver without polluting their code with SAMRAI details.
 *
 *
 * The MultiPhysicsIntegrator will use abstract ITransaction objects mainly to:
 *  - regrid a level (see regrid method), without knowing which data has to be set on the reggrided
 * level and from which data it is obtained.
 *  - intiallize a level (see initLevel method) without knowing which data is to be initialized
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
     * This routine will be typically called by the MultiPhysicsIntegrator while in
     * initializeLevelData. Typical internal resources owned by an ITransactions are temporal
     * temporary versions of model variable used for time interpolations. Other kind of internal
     * resources can be temporary variables used in multi-model transactions involving quantities of
     * different nature
     *
     * @param patch is the samrai Patch onto which we want to allocate the data
     * @param allocateTime is the time at which data is allocated
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const = 0;


    /**
     * @brief setup is used to tell the transaction which data will be communicated
     * with coarser and finer levels. This method is called only once when the Transaction
     * is created. Only information regarding which data will be communicated is given, no
     * information about which levels are concerned because that one needs to be updated
     * regularly as the patch hierarchy evolves.
     *
     * setup needs information from both the coarse and fine levels involved in the Transaction
     * as these levels may involve different models and solvers. The method takes ownership of
     * abstract ITransactionInfo structures. Concrete implementation of setup, e.g.
     * HybridTransaction will need to cast these into their derived class of ITransactionInfo, e.g.
     * HybridTransactionInfo.
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
    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       std::unique_ptr<ITransactionInfo> fromFinerInfo)
        = 0;



    /**
     * @brief setLevel tells the transaction which level in the hierarchy data will be communicated
     * to. When the transaction involves communication from coarser grid, it is from levelNumber-1.
     * This method needs to be called whenever initializing a new level in the hierarchy so that a
     * valid Transaction be used at this level. The ITransaction need to be setup before The
     * levelnumber needs to be the one of a level in the given hierarchy.
     */
    virtual void setLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
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


    virtual void initLevel(int const iLevel, double const initDataTime) const = 0;


    virtual void firstStep(PhysicalModel const& model) = 0;


    virtual void lastStep(PhysicalModel const& model) = 0;



    virtual std::string fineModelName() const = 0;

    virtual std::string coarseModelName() const = 0;



    virtual ~ITransaction() = default;
};



bool areCompatible(ITransaction const& transaction, PhysicalModel const& model);



} // namespace PHARE
#endif
