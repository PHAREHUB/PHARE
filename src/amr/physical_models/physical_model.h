


#ifndef PHARE_PHYSICAL_MODEL_H
#define PHARE_PHYSICAL_MODEL_H

#include <memory>
#include <string>

#include <SAMRAI/hier/Patch.h>

#include "evolution/messengers/messenger_info.h"

namespace PHARE
{
/**
 * @brief The PhysicalState class is an interface for concrete states holding data associated with a
 * concrete IPhysicalModel
 */
class IPhysicalState
{
};



/**
 * @brief The IPhysicalModel class represents an interface for manipulating physical quantities
 * governed by a set of equations on an AMR patch hierarchy.
 */
class IPhysicalModel
{
protected:
    IPhysicalModel(std::string modelName)
        : name_{std::move(modelName)}
    {
    }


    std::string name_;

public:
    std::string name() const { return name_; }



    /**
     * @brief allocate must be implemented by concrete subclasses to allocate the model quantities
     * onto the given patch at the given allocateTime.
     *
     * Generally, a concrete IPhysicalModel will hold a concrete PhysicalState and a
     * ResourcesManager, and this method will call the ResourcesManager allocate() method passing it
     * the quantities in the PhysicalState.
     *
     * This method is typically called in the MultiPhysicsIntegrator when initializing level data
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) = 0;



    /**
     * @brief fillMessengerInfo mut be implemented by concrete subclasses. The method is called by
     * the MessengerRegistration class to register the quantities the model needs to be initialized,
     * filled on ghosts.
     */
    virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const = 0;




    virtual ~IPhysicalModel() = default;
};



} // namespace PHARE

#endif
