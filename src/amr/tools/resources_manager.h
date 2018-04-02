#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_H
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_H

#include "hybrid/hybrid_quantities.h"
#include "resources_guards.h"
#include "resources_manager_utilities.h"

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/VariableDatabase.h>


#include <map>


namespace PHARE
{
/**
 * \brief FieldType gathers field type traits for ResourceManager to know
 * how to work with Field objects
 **/
template<typename Impl>
struct FieldType
{
    using field_impl = typename Impl::field_impl;
};



/**
 * \brief Similarly to FieldType, ParticleType gathers particle type traits for
 * the ResourceManager to work with Particle objects.
 */
template<typename Impl>
struct ParticlesType
{
    using particles_impl = typename Impl::particles_impl;
};


/**
 * \brief ResourcesInfo gathers the information associated to a resource
 * that is necessary for the ResourceManager to retrieve data.
 */
struct ResourcesInfo
{
    std::shared_ptr<SAMRAI::hier::Variable> variable;
    int id;
};

/** \brief ResourcesManager is an adapter between PHARE objects that manipulate
 * data on patches, and the SAMRAI variable database system, storing the data.
 * It is used by PHARE to register data to the samrai system and to get access to it
 * whenever needed, without having to know the SAMRAI database system.
 *
 * Objects registering and retrieving data through the ResourcesManager are called
 * ResourcesUser. A ResourcesUser needs to satisfy a specific interface to work with
 * the ResourcesManager.
 *
 * There are only two kinds of Resources that can be registered to the ResourcesManager
 *
 * - resources defined as a Field
 * - resources defined as a ParticleArray
 *
 * Several kinds of ResourcesUser can register their resources to the ResourcesManager
 * and are identified by type traits.
 *
 * - those having Field resources only (trait isFieldType)
 * - those having ParticleArray resources only (trait isParticleType)
 * - those having subResources, i.e. objects that do not hold Resources themselves
 *   but rather hold objects that hold resources. (trait isSubResourcesType)
 *
 * and any combination of those, for instance the following traits are used:
 *  - isParticleAndFieldType
 *  - isFieldAndSubResourcesType
 *
 * the following traits are not used yet since no data structure in PHARE requires it:
 * - isParticlesAndSubResourcesType
 * - isAllType (Field, Particle and SubResources)
 *
 *
 * for exemple:
 *
 * - VecField is a ResourcesUser with the isFieldType trait because it only has
 * Field resources.
 * - Electromag is a ResourcesUser with a isSubResourcesType trait because it
 * does not hold any Field or ParticleArray, but holds 2 VecFields that hold Field objects.
 *
 * ResourcesManager basically offers 4 kinds of methods in its interface. Three
 * methods aim to be used by code manipulating ResourceUsers, typically the
 * code that manipulate data and algorithms:
 *
 * - registerResources() : these methods are used by code that wants to register
 * resources of a ResourcesUser.
 *
 * - setResources() : these methods are used by code that wants to retrieve a pointer
 * to an already registered resource. After calling this method passing a ResourcesUser
 * and a Patch, the ResourcesUser has its internal pointers pointing to the resource
 * lying on the given Patch. These methods are typically used while looping over the
 * patches of a PatchLevel, to set the resources of a ResourcesUser to those lying on
 * the current patch. One issue is that ResourcesUser internal pointers
 * will point to the data of the last Patch it has been set to. This is why this
 * method is typically not used directly. Rather, it is better to use a ResourcesGuard
 * that will set the resources of ResourcesUsers in scope and reset their pointers
 * upon leaving the scope.
 *
 * - createResourcesGuard() : used to create a ResourcesGuard.
 *
 * One method is rather to be used by SAMRAI derived classes:
 *
 * - allocate()
 *
 */
class ResourcesManager
{
public:
    ResourcesManager() = delete;

    ResourcesManager(std::string const& layoutType, SAMRAI::tbox::Dimension const& dim)
        : gridLayout_{layoutType}
        , variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , context_{variableDatabase_->getContext(contextName_)}
        , dimension_{dim}
    {
    }


    ResourcesManager(ResourcesManager const&) = delete;
    ResourcesManager(ResourcesManager&&)      = delete;
    ResourcesManager& operator=(ResourcesManager const& source) = delete;
    ResourcesManager& operator=(ResourcesManager&&) = delete;


    /**
     * \brief register the resources of a ResourcesUser that only has Field resources
     *
     * isFieldType<ResourcesUser> SFINAEs out this specialization if ResourcesUser
     * does not have getFieldNamesAndQuantities(). The same kind of compile time
     * introspection is used for all specializations.
     *
     * \param[in] obj
     */
    template<typename ResourcesUser, isFieldType<ResourcesUser> = dummy::value>
    void registerResources(ResourcesUser const& obj)
    {
        registerResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities());
    }




    template<typename ResourcesUser, isParticlesType<ResourcesUser> = dummy::value>
    void registerResources(ResourcesUser const& obj)
    {
        registerResources_(obj, ParticlesType<ResourcesUser>{},
                           obj.getParticlesNamesAndQuantities());
    }




    template<typename ResourcesUser, isSubResourcesType<ResourcesUser> = dummy::value>
    void registerResources(ResourcesUser const& obj)
    {
        registerSubResources_(obj.getSubResourcesObject());
    }




    template<typename ResourcesUser, isFieldAndParticlesType<ResourcesUser> = dummy::value>
    void registerResources(ResourcesUser const& obj)
    {
        registerResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities());
        registerResources_(obj, ParticlesType<ResourcesUser>{},
                           obj.getParticlesNamesAndQuantities());
    }




    template<typename ResourcesUser, isFieldAndSubResourcesType<ResourcesUser> = dummy::value>
    void registerResources(ResourcesUser const& obj)
    {
        registerResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities());
        registerSubResources_(obj.getSubResourcesObject());
    }




#if PARTICLES_AND_SUB_RESOURCES_EXIST
    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj,
                           isParticlesAndSubResourcesType<ResourcesUser> = std::true_type{})
    {
    }
#endif
#if FIELD_PARTICLES_AND_SUB_RESOURCES_EXIST
    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj, isAllType<ResourcesUser> = std::true_type{})
    {
    }
#endif



    template<typename ResourcesUser, typename NullOrResourcePtr,
             isFieldType<ResourcesUser> = dummy::value>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch) const
    {
        setResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities(), patch,
                      nullOrResourcePtr);
    }



    template<typename ResourcesUser, typename NullOrResourcePtr,
             isParticlesType<ResourcesUser> = dummy::value>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch) const
    {
        setResources_(obj, ParticlesType<ResourcesUser>{}, obj.getParticleArrayNames(), patch,
                      nullOrResourcePtr);
    }


    template<typename ResourcesUser, typename NullOrResourcePtr,
             isSubResourcesType<ResourcesUser> = dummy::value>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch) const
    {
        setResourcesVariadic_(nullOrResourcePtr, patch, obj.getSubResourcesObject());
    }



    template<typename ResourcesUser, typename NullOrResourcePtr,
             isFieldAndParticlesType<ResourcesUser> = dummy::value>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch) const
    {
        setResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities(), patch,
                      nullOrResourcePtr);

        setResources_(obj, ParticlesType<ResourcesUser>{}, obj.getParticleArrayNames(), patch,
                      nullOrResourcePtr);
    }



    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isFieldAndSubResourcesType<ResourcesUser> = std::true_type{}) const
    {
        setResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities(), patch,
                      nullOrResourcePtr);

        setResourcesVariadic_(nullOrResourcePtr, patch, obj.getSubResourcesObject());
    }


#if PARTICLES_AND_SUB_RESOURCES_EXIST
    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isParticlesAndSubResourcesType<ResourcesUser> = std::true_type{}) const
    {
    }
#endif
#if FIELD_PARTICLES_AND_SUB_RESOURCES_EXIST
    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isAllType<ResourcesUser> = std::true_type{}) const
    {
    }
#endif


    template<typename ResourcesUser, isFieldType<ResourcesUser> = dummy::value>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch) const
    {
        allocate_(obj, obj.getFieldNamesAndQuantities(), patch);
    }



    template<typename ResourcesUser, isParticlesType<ResourcesUser> = dummy::value>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch) const
    {
        allocate_(obj, obj.getParticleArrayNames(), patch);
    }



    template<typename ResourcesUser, isSubResourcesType<ResourcesUser> = dummy::value>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch) const
    {
        allocateVariadic_(patch, obj.getSubResourcesObject());
    }



    template<typename ResourcesUser, isFieldAndParticlesType<ResourcesUser> = dummy::value>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch) const
    {
        allocate_(obj, obj.getFieldNamesAndQuantities(), patch);
        allocate_(obj, obj.getParticleArrayNames(), patch);
    }


    template<typename ResourcesUser, isFieldAndSubResourcesType<ResourcesUser> = dummy::value>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch) const
    {
        allocate_(obj, obj.getFieldNamesAndQuantities(), patch);
        allocateVariadic_(patch, obj.getSubResourcesObject());
    }




#if PARTICLES_AND_SUB_RESOURCES_EXIST
    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isParticlesAndSubResourcesType<ResourcesUser> = std::true_type{}) const
    {
    }
#endif
#if FIELD_PARTICLES_AND_SUB_RESOURCES_EXIST
    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isAllType<ResourcesUser> = std::true_type{}) const
    {
    }
#endif


    // This is a helper that allow use simple as
    // auto guard = createResourcesGuards(patch, obj1, obj2, ...);

    template<typename... ResourcesUsers>
    constexpr ResourcesGuards<ResourcesManager, ResourcesUsers...>
    createResourcesGuards(SAMRAI::hier::Patch const& patch, ResourcesUsers&... resoucesUsers)
    {
        return ResourcesGuards<ResourcesManager, ResourcesUsers...>{patch, *this, resoucesUsers...};
    }

private:
    // The function getPointer_ is the one that depending on NullOrResourcePtr will choose to return
    // the real pointer or a nullptr to the correct type. Note that it is mandatory to return
    // a nullptr of correct type, otherwise the call of ResourcesUser.setResources(name,nullptr)
    // will be undefined. To accomplish that, ResourcesType will have to provides the type of the
    // internal data (which is known by core)

    /** \brief Return a pointer to the patch data instantied in the patch
     */
    template<typename ResourceType>
    auto getResourcePointer_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                             SAMRAI::hier::Patch const& patch) const
    {
        auto patchData = patch.getPatchData(resourcesVariableInfo.variable, context_);
        return (std::dynamic_pointer_cast<typename ResourceType::data_type>(patchData))
            ->getPointer();
    }




    template<typename ResourceType>
    constexpr auto static getNullPointer_(ResourceType resourceType)
    {
        return static_cast<typename ResourceType::internal_type_ptr>(nullptr);
    }



    template<typename ResourceType, typename NullOrResourcePtr,
             ifWantRessourcePointer<NullOrResourcePtr> = dummy::value>
    auto getPointer_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                     SAMRAI::hier::Patch const& patch, NullOrResourcePtr which) const
    {
        return getResourcePointer_(resourceType, resourcesVariableInfo, patch);
    }




    template<typename ResourceType, typename NullOrResourcePtr,
             ifWantNullPointer<NullOrResourcePtr> = dummy::value>
    auto getPointer_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                     SAMRAI::hier::Patch const& patch, NullOrResourcePtr which) const
    {
        (void)resourcesVariableInfo;
        (void)patch;
        return getNullPointer_(resourceType);
    }



    template<typename ResourcesUser, typename... ResourcesUserTails>
    void registerSubResources_(ResourcesUser& resource, ResourcesUserTails&... resources)
    {
        registerResources(resource);
        registerSubResources_(resources...);
    }




    template<typename NullOrResourcePtr, typename ResourcesUser, typename... ResourcesUserTails>
    void setResourcesVariadic_(NullOrResourcePtr nullOrResourcePtr,
                               SAMRAI::hier::Patch const& patch, ResourcesUser obj,
                               ResourcesUserTails&... objs) const
    {
        setResources(obj, nullOrResourcePtr, patch);
        setResourcesVariadic_(nullOrResourcePtr, patch, objs...);
    }




    template<typename ResourcesUser, typename... ResourcesUserTails>
    void allocateVariadic_(SAMRAI::hier::Patch const& patch, ResourcesUser obj,
                           ResourcesUserTails&... objs)
    {
        allocate(obj, patch);
        allocateVariadic_(patch, objs...);
    }




    /** \brief Register a collection of resources
     */
    template<typename ResourcesUser, typename ResourcesType>
    void registerResources_(ResourcesUser const& obj, ResourcesType resourceType,
                            typename ResourcesUser::names_and_quantities const& namesAndQuantities)
    {
        for (auto const& nameAndQuantity : namesAndQuantities)
        {
            std::string const& resourcesName = nameAndQuantity.first;
            auto const& qty                  = nameAndQuantity.second;

            ResourcesInfo resources;
            resources.variable = std::make_shared<typename ResourcesType::variable_type>(
                dimension_, resourcesName, qty, gridLayout_);

            resources.id = variableDatabase_->registerVariableAndContext(
                resources.variable, context_, SAMRAI::hier::IntVector::getZero(dimension_));

            nameToressourceInfo_[resourcesName] = resources;
        }
    }



    /** \brief setResources_ aims at setting ResourcesUser pointers to the
     * appropriate data on the Patch or to reset them to nullptr.
     */
    template<typename ResourcesUser, typename ResourcesType, typename NullOrResourcePtr>
    void setResources_(ResourcesUser& obj, ResourcesType resourceType,
                       typename ResourcesUser::names_and_quantities const& namesAndQuantities,
                       SAMRAI::hier::Patch const& patch, NullOrResourcePtr nullOrResourcePtr) const
    {
        for (auto const& nameAndQuantity : namesAndQuantities)
        {
            std::string const& resourcesName = nameAndQuantity.first;
            auto const& resourceInfoIt       = nameToressourceInfo_.find(resourcesName);
            if (resourceInfoIt != nameToressourceInfo_.end())
            {
                auto data
                    = getPointer_(resourceType, resourceInfoIt->second, patch, nullOrResourcePtr);
                obj.setResources(resourcesName, data);
            }
            else
            {
                throw std::runtime_error("Resources not found !");
            }
        }
    }




    //! \brief Allocate the data on the given level
    template<typename ResourcesUser>
    void allocate_(ResourcesUser const& obj,
                   typename ResourcesUser::names_and_quantities const& namesAndQuantities,
                   SAMRAI::hier::Patch& patch) const
    {
        for (auto const& nameAndQuantity : namesAndQuantities)
        {
            std::string const& resourcesName  = nameAndQuantity.first;
            auto const& resourceVariablesInfo = nameToressourceInfo_.find(resourcesName);
            if (resourceVariablesInfo != nameToressourceInfo_.end())
            {
                patch.allocatePatchData(resourceVariablesInfo->second.id);
            }
            else
            {
                throw std::runtime_error("Resources not found !");
            }
        }
    }

    std::string contextName_{"default"};
    std::string gridLayout_{"yee"};
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    SAMRAI::tbox::Dimension dimension_;
    std::map<std::string, ResourcesInfo> nameToressourceInfo_;
}; // namespace PHARE

} // namespace PHARE
#endif
