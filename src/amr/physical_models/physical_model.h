


#ifndef PHARE_PHYSICAL_MODEL_H
#define PHARE_PHYSICAL_MODEL_H

#include <memory>
#include <string>

#include <SAMRAI/hier/Patch.h>

#include "evolution/transactions/transaction_info.h"

namespace PHARE
{
class PhysicalState
{
};


class PhysicalModel
{
protected:
    PhysicalModel(std::string modelName)
        : name_{std::move(modelName)}
    {
    }


    std::string name_;

public:
    virtual std::unique_ptr<ITransactionInfo> transactionInfo() const = 0;
    std::string name() const { return name_; }
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) = 0;

    virtual void fillTransactionInfo(std::unique_ptr<ITransactionInfo> const& info) const = 0;

    virtual ~PhysicalModel() = default;
};



} // namespace PHARE

#endif
