
#ifndef PHARE_SOLVER_HPP
#define PHARE_SOLVER_HPP

#include <string>

#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/physical_models/physical_model.hpp"



namespace PHARE
{
namespace solver
{
    /**
     * @brief The ISolver is an interface for a solver used by the MultiPHysicsIntegrator.
     *
     * The main interest of this class is to provide the MultiPhysicsIntegrator with the method
     * advanceLevel().
     *
     */
    template<typename AMR_Types>
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
         * @brief registerResources is used to register the solver quantities that need to be
         * defined on a Patch. The quantities are registered to the ResourcesManager of the given
         * IPhysicalModel
         * @param model
         */
        virtual void registerResources(IPhysicalModel<AMR_Types>& model) = 0;




        /**
         * @brief fillMessengerInfo fills the IMessengerInfo with the names of the ISolver
         * quantities that need to be communicated by a IMessenger.
         * @param info
         */
        virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const = 0;




        /**
         * @brief advanceLevel advances the given level from t to t+dt
         */
        virtual void advanceLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                  int const levelNumber, IPhysicalModel<AMR_Types>& model,
                                  amr::IMessenger<IPhysicalModel<AMR_Types>>& fromCoarser,
                                  const double currentTime, const double newTime)
            = 0;




        /**
         * @brief allocate is used to allocate ISolver variables previously registered to the
         * ResourcesManager of the given model, onto the given Patch, at the given time.
         */
        virtual void allocate(IPhysicalModel<AMR_Types>& model, SAMRAI::hier::Patch& patch,
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


    template<typename AMR_Types>
    bool areCompatible(amr::IMessenger<IPhysicalModel<AMR_Types>> const& messenger,
                       ISolver<AMR_Types> const& solver)
    {
        return solver.modelName() == messenger.fineModelName()
               || solver.modelName() == messenger.coarseModelName();
    }


    /**
     * @brief areCompatible returns true if the model name is equal to the solver modelname
     */
    template<typename AMR_Types>
    bool areCompatible(IPhysicalModel<AMR_Types> const& model, ISolver<AMR_Types> const& solver)
    {
        return model.name() == solver.modelName();
    }


} // namespace solver
} // namespace PHARE

#endif
