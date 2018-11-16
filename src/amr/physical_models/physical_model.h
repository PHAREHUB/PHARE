


#ifndef PHARE_PHYSICAL_MODEL_H
#define PHARE_PHYSICAL_MODEL_H

#include <memory>
#include <string>

#include <SAMRAI/hier/Patch.h>

#include "evolution/messengers/messenger_info.h"

namespace PHARE
{
class PhysicalState
{
};


class IPhysicalModel
{
protected:
    IPhysicalModel(std::string modelName)
        : name_{std::move(modelName)}
    {
    }


    std::string name_;

public:
    virtual std::unique_ptr<IMessengerInfo> messengerInfo() const = 0;
    std::string name() const { return name_; }
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) = 0;

    virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const = 0;

    virtual ~IPhysicalModel() = default;
};



} // namespace PHARE

#endif
