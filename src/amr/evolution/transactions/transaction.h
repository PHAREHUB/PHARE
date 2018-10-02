
#ifndef PHARE_TRANSACTION_H
#define PHARE_TRANSACTION_H

#include <string>

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include "physical_models/physical_model.h"
#include "transaction_info.h"


namespace PHARE
{
class ITransaction
{
public:
    virtual std::string name() = 0;

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() = 0;
    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner()   = 0;

    virtual void allocate(PhysicalModel const& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const = 0;

    virtual void initialize(PhysicalModel const& destModel, PhysicalModel const& srcModel) = 0;

    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       std::unique_ptr<ITransactionInfo> fromFinerInfo)
        = 0;

    virtual void update(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        int const levelNumber)
        = 0;


    virtual void
    regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy, const int levelNumber,
           std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel, double const initDataTime)
        = 0;


    virtual void initLevel(int const iLevel, double const initDataTime) const = 0;

    virtual std::string fineModelName() const = 0;

    virtual std::string coarseModelName() const = 0;



    virtual ~ITransaction() = default;
};



bool areCompatible(ITransaction const& transaction, PhysicalModel const& model);



} // namespace PHARE
#endif
