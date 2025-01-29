#ifndef PHARE_PHYSICAL_MODEL_HPP
#define PHARE_PHYSICAL_MODEL_HPP

#include <memory>
#include <string>

#include "amr/messengers/messenger_info.hpp"
#include "amr/resources_manager/resources_manager.hpp"

namespace PHARE
{
namespace solver
{
    /**
     * @brief The IPhysicalModel class represents an interface for manipulating physical quantities
     * governed by a set of equations on an AMR patch hierarchy.
     */
    template<typename AMR_Types>
    class IPhysicalModel
    {
    protected:
        IPhysicalModel(std::string modelName)
            : name_{std::move(modelName)}
        {
        }


        std::string name_;

    public:
        using patch_t = typename AMR_Types::patch_t;
        using level_t = typename AMR_Types::level_t;

        std::string name() const { return name_; }



        /**
         * @brief initialize is used to initialize data
         */
        virtual void initialize(level_t& level) = 0;


        virtual amr::IResourcesManager& resources_manager() const = 0;


        /**
         * @brief allocate must be implemented by concrete subclasses to allocate the model
         * quantities onto the given patch at the given allocateTime.
         *
         * Generally, a concrete IPhysicalModel will hold a concrete PhysicalState and a
         * ResourcesManager, and this method will call the ResourcesManager allocate() method
         * passing it the quantities in the PhysicalState.
         *
         * This method is typically called in the MultiPhysicsIntegrator when initializing level
         * data
         */
        virtual void allocate(patch_t& patch, double const allocateTime) = 0;



        /**
         * @brief fillMessengerInfo mut be implemented by concrete subclasses. The method is called
         * by the MessengerRegistration class to register the quantities the model needs to be
         * initialized, filled on ghosts.
         */
        virtual void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const = 0;




        virtual ~IPhysicalModel() = default;
    };
} // namespace solver


} // namespace PHARE

#endif
