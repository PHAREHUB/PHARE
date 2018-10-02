
#ifndef PHARE_SOLVER_H
#define PHARE_SOLVER_H

#include <string>

#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>

#include "evolution/transactions/transaction.h"
#include "evolution/transactions/transaction_info.h"
#include "physical_models/physical_model.h"



namespace PHARE
{
class ISolver
{
public:
    std::string name() const { return solverName; }

    virtual std::string modelName() const = 0;

    virtual void registerResources(PhysicalModel& model) = 0;

    virtual void fillTransactionInfo(std::unique_ptr<ITransactionInfo> const& info) const = 0;


    virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              int const levelNumber, PhysicalModel& model,
                              ITransaction& fromCoarser, const double currentTime,
                              const double newTime)
        = 0;


    virtual void allocate(PhysicalModel& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const = 0;

    virtual ~ISolver() = default;


protected:
    explicit ISolver(std::string name)
        : solverName{std::move(name)}
    {
    }
    std::string solverName;
};



bool areCompatible(ITransaction const& transaction, ISolver const& solver);

bool areCompatible(PhysicalModel const& model, ISolver const& solver);



} // namespace PHARE

#endif
