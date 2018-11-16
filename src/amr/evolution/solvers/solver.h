
#ifndef PHARE_SOLVER_H
#define PHARE_SOLVER_H

#include <string>

#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>

#include "evolution/messengers/messenger.h"
#include "evolution/messengers/messenger_info.h"
#include "physical_models/physical_model.h"



namespace PHARE
{
class ISolver
{
public:
    std::string name() const { return solverName; }

    virtual std::string modelName() const = 0;

    virtual void registerResources(IPhysicalModel& model) = 0;

    virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const = 0;


    virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              int const levelNumber, IPhysicalModel& model, IMessenger& fromCoarser,
                              const double currentTime, const double newTime)
        = 0;


    virtual void allocate(IPhysicalModel& model, SAMRAI::hier::Patch& patch,
                          double const allocateTime) const = 0;

    virtual ~ISolver() = default;


protected:
    explicit ISolver(std::string name)
        : solverName{std::move(name)}
    {
    }
    std::string solverName;
};



bool areCompatible(IMessenger const& messenger, ISolver const& solver);

bool areCompatible(IPhysicalModel const& model, ISolver const& solver);



} // namespace PHARE

#endif
