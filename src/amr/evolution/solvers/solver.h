
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
/**
 * @brief The ISolver is an interface for a solver used by the MultiPHysicsIntegrator.
 *
 * The main interest of this class is to provide the MultiPhysicsIntegrator with the method
 * advanceLevel().
 *
 */
class ISolver
{
public:
    /**
     * @brief return the name of the ISolver
     */
    std::string name() const { return solverName; }



    /**
     * @brief return the name of the model the ISolver is compatible with
     */
    virtual std::string modelName() const = 0;



    /**
     * @brief registerResources is used to register the solver quantities that need to be defined on
     * a Patch. The quantities are registered to the ResourcesManager of the given IPhysicalModel
     * @param model
     */
    virtual void registerResources(IPhysicalModel& model) = 0;




    /**
     * @brief fillMessengerInfo fills the IMessengerInfo with the names of the ISolver quantities
     * that need to be communicated by a IMessenger.
     * @param info
     */
    virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const = 0;




    /**
     * @brief advanceLevel advances the given level from t to t+dt
     */
    virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              int const levelNumber, IPhysicalModel& model, IMessenger& fromCoarser,
                              const double currentTime, const double newTime)
        = 0;




    /**
     * @brief allocate is used to allocate ISolver variables previously registered to the
     * ResourcesManager of the given model, onto the given Patch, at the given time.
     */
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


/**
 * @brief areCompatible returns true if the solver.modelname() equals the messenger name
 */
bool areCompatible(IMessenger const& messenger, ISolver const& solver);


/**
 * @brief areCompatible returns true if the model name is equal to the solver modelname
 */
bool areCompatible(IPhysicalModel const& model, ISolver const& solver);



} // namespace PHARE

#endif
