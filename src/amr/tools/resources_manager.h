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
 * \brief FieldType is a structure that will contain some
 * type information related to Field
 *
 **/
template<typename Impl>
struct FieldType
{
    using field_impl = typename Impl::field_impl;
};
/**
 * \brief ParticlesType is a structure that will contain some
 * type information related to Particles
 */
template<typename Impl>
struct ParticlesType
{
    using particles_impl = typename Impl::particles_impl;
};

struct ResourcesInfo
{
    std::shared_ptr<SAMRAI::hier::Variable> variable;
    int id;
};

/** \brief ResourcesManager allow object that need to manipulate data relative
 * to a patch, to register them and when needed to give access to the data
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
     * \brief Overload for ResourcesUser that only have field
     * \param[in] obj
     */
    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj, isFieldType<ResourcesUser> = std::true_type{})
    {
        registerResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities());
    }




    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj,
                           isParticlesType<ResourcesUser> = std::true_type{})
    {
        registerResources_(obj, ParticlesType<ResourcesUser>{},
                           obj.getParticlesNamesAndQuantities());
    }




    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj,
                           isSubResourcesType<ResourcesUser> = std::true_type{})
    {
        registerResourcesVariadic_(obj.getSubResourcesObject());
    }




    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj,
                           isFieldAndParticlesType<ResourcesUser> = std::true_type{})
    {
        registerResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities());
        registerResources_(obj, ParticlesType<ResourcesUser>{},
                           obj.getParticlesNamesAndQuantities());
    }




    template<typename ResourcesUser>
    void registerResources(ResourcesUser const& obj,
                           isFieldAndSubResourcesType<ResourcesUser> = std::true_type{})
    {
        registerResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities());
        registerResourcesVariadic_(obj.getSubResourcesObject());
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


    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isFieldType<ResourcesUser> = std::true_type{}) const
    {
        setResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities(), patch,
                      nullOrResourcePtr);
    }




    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isParticlesType<ResourcesUser> = std::true_type{}) const
    {
        setResources_(obj, ParticlesType<ResourcesUser>{}, obj.getParticlesdNamesAndQuantities(),
                      patch, nullOrResourcePtr);
    }




    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isSubResourcesType<ResourcesUser> = std::true_type{}) const
    {
        setResourcesVariadic_(nullOrResourcePtr, patch, obj.getSubResourcesObject());
    }



    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                      SAMRAI::hier::Patch const& patch,
                      isFieldAndParticlesType<ResourcesUser> = std::true_type{}) const
    {
        setResources_(obj, FieldType<ResourcesUser>{}, obj.getFieldNamesAndQuantities(), patch,
                      nullOrResourcePtr);

        setResources_(obj, ParticlesType<ResourcesUser>{}, obj.getParticlesdNamesAndQuantities(),
                      patch, nullOrResourcePtr);
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


    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isFieldType<ResourcesUser> = std::true_type{}) const
    {
        allocate_(obj, obj.getFieldNamesAndQuantities(), patch);
    }



    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isParticlesType<ResourcesUser> = std::true_type{}) const
    {
        allocate_(obj, obj.getParticlesdNamesAndQuantities(), patch);
    }



    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isSubResourcesType<ResourcesUser> = std::true_type{}) const
    {
        allocateVariadic_(patch, obj.getSubResourcesObject());
    }



    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isFieldAndParticlesType<ResourcesUser> = std::true_type{}) const
    {
        allocate_(obj, obj.getFieldNamesAndQuantities(), patch);
        allocate_(obj, obj.getParticlesdNamesAndQuantities(), patch);
    }


    template<typename ResourcesUser>
    void allocate(ResourcesUser const& obj, SAMRAI::hier::Patch& patch,
                  isFieldAndSubResourcesType<ResourcesUser> = std::true_type{}) const
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



    template<typename ResourceType, typename NullOrResourcePtr>
    auto getPointer_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                     SAMRAI::hier::Patch const& patch, NullOrResourcePtr which,
                     ifWantRessourcePointer<NullOrResourcePtr> = std::true_type{}) const
    {
        return getResourcePointer_(resourceType, resourcesVariableInfo, patch);
    }




    template<typename ResourceType, typename NullOrResourcePtr>
    auto getPointer_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                     SAMRAI::hier::Patch const& patch, NullOrResourcePtr which,
                     ifWantNullPointer<NullOrResourcePtr> = std::true_type{}) const
    {
        (void)resourcesVariableInfo;
        (void)patch;
        return getNullPointer_(resourceType);
    }



    template<typename ResourcesUser, typename... ResourcesUserTails>
    void registerResourcesVariadic_(ResourcesUser& resource, ResourcesUserTails&... resources)
    {
        registerResources(resource);
        registerResourcesVariadic_(resources...);
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




    /** \brief Set the resources of an object
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
};

} // namespace PHARE
#endif
