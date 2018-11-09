#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_H
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_H

#include "field_resource.h"
#include "hybrid/hybrid_quantities.h"
#include "particle_resource.h"
#include "resources_guards.h"
#include "resources_manager_utilities.h"


#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/VariableDatabase.h>


#include <map>
#include <optional>


namespace PHARE
{
/**
 * \brief FieldType gathers field type traits for ResourceManager to know
 * how to work with Field objects
 **/




template<typename GridLayoutT, typename ResourcesUser, typename ResourcesType>
using isUserFieldType = std::is_same<std::remove_reference_t<ResourcesType>,
                                     UserFieldType<GridLayoutT, ResourcesUser>>;




template<typename ResourcesUser, typename ResourcesType>
using isUserParticleType = std::is_same<typename std::remove_reference<ResourcesType>::type,
                                        UserParticleType<ResourcesUser>>;




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
 * - those having Field resources only (trait has_field)
 * - those having ParticleArray resources only (trait has_particle)
 * - those having subResources, i.e. objects that do not hold Resources themselves
 *   but rather hold objects that hold resources. (trait has_sub_resources)
 *
 *
 * for exemple:
 *
 * - VecField is a ResourcesUser for which the has_field trait is true
 * - Electromag is a ResourcesUser for which the has_sub_resources trait  is true because it
 *  holds 2 VecFields that hold Field objects.
 *
 * ResourcesManager basically offers 4 kinds of methods in its interface. Three
 * methods aim to be used by code manipulating ResourceUsers, typically the
 * code that manipulate data and algorithms:
 *
 * - registerResources() : these methods are used by code that wants to register
 * resources of a ResourcesUser.
 *
 * - setResources_() : these methods are used by code that wants to retrieve a pointer
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
template<typename GridLayoutT>
class ResourcesManager
{
public:
    ResourcesManager()
        : variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , context_{variableDatabase_->getContext(contextName_)}
        , dimension_{SAMRAI::tbox::Dimension{GridLayoutT::dimension}}
    {
    }


    ResourcesManager(ResourcesManager const&) = delete;
    ResourcesManager(ResourcesManager&&)      = delete;
    ResourcesManager& operator=(ResourcesManager const& source) = delete;
    ResourcesManager& operator=(ResourcesManager&&) = delete;


    /** @brief registerResources takes a ResourcesUser to register its resources
     * to the SAMRAI data system.
     *
     * the function asks the ResourcesUser informations (at least the name) of the
     * buffers to register. This information is used to know which kind of
     * SAMRAI patchdata & variable to create, but is also stored to retrieve
     * this data when the ResourcesUser is later going to ask for it (via the BufferGuard)
     *
     * The function registers all FieldData for ResourcesUser that have Fields, all ParticleData for
     * ResourcesUser that have particle Arrays. In case the ResourcesUser has sub-resources
     * we ask for them in a tuple, and recursively call registerResources() for all of the unpacked
     * elements
     */
    template<typename ResourcesUser>
    void registerResources(ResourcesUser& obj)
    {
        if constexpr (has_field<ResourcesUser>::value)
        {
            registerResources_<ResourcesUser, UserFieldType<GridLayoutT, ResourcesUser>>(obj);
        }

        if constexpr (has_particles<ResourcesUser>::value)
        {
            registerResources_<ResourcesUser, UserParticleType<ResourcesUser>>(obj);
        }


        if constexpr (has_runtime_subresourceuser_list<ResourcesUser>::value)
        {
            auto&& resourcesUsers = obj.getRunTimeResourcesUserList();
            for (auto& resourcesUser : resourcesUsers)
            {
                this->registerResources(resourcesUser);
            }
        }


        if constexpr (has_compiletime_subresourcesuser_list<ResourcesUser>::value)
        {
            // get a tuple here
            auto&& subResources = obj.getCompileTimeResourcesUserList();

            // unpack the tuple subResources and apply for each element registerResources()
            // (recursively)
            std::apply(
                [this](auto&... subResource) { (this->registerResources(subResource), ...); },
                subResources);
        }
    }



    /** @brief allocate the appropriate PatchDatas on the Patch for the ResourcesUser
     *
     * The function allocates all FieldData for ResourcesUser that have Fields, all ParticleData for
     * ResourcesUser that have particle Arrays. In case the ResourcesUser has sub-resources
     * we ask for them in a tuple, and recursively call allocate() for all of the unpacked elements
     */
    template<typename ResourcesUser>
    void allocate(ResourcesUser& obj, SAMRAI::hier::Patch& patch, double const allocateTime) const
    {
        if constexpr (has_field<ResourcesUser>::value)
        {
            allocate_(obj, obj.getFieldNamesAndQuantities(), patch, allocateTime);
        }

        if constexpr (has_particles<ResourcesUser>::value)
        {
            allocate_(obj, obj.getParticleArrayNames(), patch, allocateTime);
        }

        if constexpr (has_runtime_subresourceuser_list<ResourcesUser>::value)
        {
            auto&& resourcesUsers = obj.getRunTimeResourcesUserList();
            for (auto& resourcesUser : resourcesUsers)
            {
                this->allocate(resourcesUser, patch, allocateTime);
            }
        }

        if constexpr (has_compiletime_subresourcesuser_list<ResourcesUser>::value)
        {
            // get a tuple here
            auto&& subResources = obj.getCompileTimeResourcesUserList();

            // unpack the tuple subResources and apply for each element registerResources()
            std::apply([this, &patch, allocateTime](auto&... subResource) //
                       { (this->allocate(subResource, patch, allocateTime), ...); },
                       subResources);
        }
    }




    /** \brief make a ResourceGuard to set resources to all passed objects
     *
     * This is a helper that allow use simple as
     * auto guard = createResourcesGuards(patch, obj1, obj2, ...);
     */
    template<typename... ResourcesUsers>
    constexpr ResourcesGuard<ResourcesManager, ResourcesUsers...>
    makeResourcesGuard(SAMRAI::hier::Patch const& patch, ResourcesUsers&... resourcesUsers)
    {
        return ResourcesGuard<ResourcesManager, ResourcesUsers...>{patch, *this, resourcesUsers...};
    }


    template<typename ResourcesUser>
    void setTime(ResourcesUser& obj, SAMRAI::hier::Patch const& patch, double time) const
    {
        auto IDs = getIDs(obj);
        for (auto const& id : IDs)
        {
            auto patchdata = patch.getPatchData(id);
            patchdata->setTime(time);
        }
    }



    /** \brief Get all the names and resources id that the resource user
     *  have registered via the ResourcesManager
     */
    template<typename ResourcesUser>
    std::vector<int> getIDs(ResourcesUser& obj) const
    {
        std::vector<int> IDs;
        this->getIDs_(obj, IDs);
        return IDs;
    }



    /** \brief Get all the names and resources id that the resource user
     *  have registered via the ResourcesManager
     */
    std::optional<int> getID(std::string const& resourceName) const
    {
        auto id = nameToResourceInfo_.find(resourceName);
        if (id != std::end(nameToResourceInfo_))
            return std::optional<int>{id->second.id};
        else
            return std::nullopt;
    }




private:
    template<typename ResourcesUser>
    void getIDs_(ResourcesUser& obj, std::vector<int>& IDs) const
    {
        if constexpr (has_field<ResourcesUser>::value)
        {
            for (auto const& properties : obj.getFieldNamesAndQuantities())
            {
                auto foundIt = nameToResourceInfo_.find(properties.name);
                if (foundIt != nameToResourceInfo_.end())
                {
                    IDs.push_back(foundIt->second.id);
                }
                else
                {
                    throw std::runtime_error("Cannot find " + properties.name);
                }
            }
        }

        if constexpr (has_particles<ResourcesUser>::value)
        {
            for (auto const& properties : obj.getParticleArrayNames())
            {
                auto foundIt = nameToResourceInfo_.find(properties.name);
                if (foundIt != nameToResourceInfo_.end())
                {
                    IDs.push_back(foundIt->second.id);
                }
                else
                {
                    throw std::runtime_error("Cannot find " + properties.name);
                }
            }
        }

        if constexpr (has_runtime_subresourceuser_list<ResourcesUser>::value)
        {
            auto&& resourcesUsers = obj.getRunTimeResourcesUserList();
            for (auto& resourcesUser : resourcesUsers)
            {
                //
                this->getIDs_(resourcesUser, IDs);
            }
        }

        if constexpr (has_compiletime_subresourcesuser_list<ResourcesUser>::value)
        {
            // get a tuple here
            auto&& subResources = obj.getCompileTimeResourcesUserList();

            // unpack the tuple subResources and apply for each element registerResources()
            std::apply(
                [this, &IDs](auto&... subResource) { (this->getIDs_(subResource, IDs), ...); },
                subResources);
        }
    }


    // The function getResourcesPointer_ is the one that depending
    // on NullOrResourcePtr will choose to return
    // the real pointer or a nullptr to the correct type.

    /** \brief Returns a pointer to the patch data instantiated in the patch
     * is used by getResourcesPointer_ when user code wants the pointer to the data
     */
    template<typename ResourceType>
    auto getPatchData_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                       SAMRAI::hier::Patch const& patch) const
    {
        auto patchData = patch.getPatchData(resourcesVariableInfo.variable, context_);
        return (std::dynamic_pointer_cast<typename ResourceType::patch_data_type>(patchData))
            ->getPointer();
    }




    /** \brief this specialization of getResourcesPointer_ is used when
     * the client code wants to get a pointer to a patch data resource
     */
    template<typename ResourceType, typename RequestedPtr>
    auto getResourcesPointer_(ResourceType resourceType, ResourcesInfo const& resourcesVariableInfo,
                              SAMRAI::hier::Patch const& patch) const
    {
        if constexpr (std::is_same_v<RequestedPtr, UseResourcePtr>)
        {
            return getPatchData_(resourceType, resourcesVariableInfo, patch);
        }

        else if constexpr (std::is_same_v<RequestedPtr, UseNullPtr>)
        {
            return static_cast<typename ResourceType::internal_type_ptr>(nullptr);
        }
    }




    template<typename ResourcesUser, typename NullOrResourcePtr>
    void setResources_(ResourcesUser& obj, NullOrResourcePtr nullOrResourcePtr,
                       SAMRAI::hier::Patch const& patch) const
    {
        if constexpr (has_field<ResourcesUser>::value)
        {
            setResourcesInternal_(obj, UserFieldType<GridLayoutT, ResourcesUser>{},
                                  obj.getFieldNamesAndQuantities(), patch, nullOrResourcePtr);
        }

        if constexpr (has_particles<ResourcesUser>::value)
        {
            setResourcesInternal_(obj, UserParticleType<ResourcesUser>{},
                                  obj.getParticleArrayNames(), patch, nullOrResourcePtr);
        }


        if constexpr (has_runtime_subresourceuser_list<ResourcesUser>::value)
        {
            auto&& resourcesUsers = obj.getRunTimeResourcesUserList();
            for (auto& resourcesUser : resourcesUsers)
            {
                this->setResources_(resourcesUser, nullOrResourcePtr, patch);
            }
        }


        if constexpr (has_compiletime_subresourcesuser_list<ResourcesUser>::value)
        {
            // get a tuple here
            auto&& subResources = obj.getCompileTimeResourcesUserList();

            // unpack the tuple subResources and apply for each element setResources_()
            std::apply(
                [this, &patch, &nullOrResourcePtr](auto&... subResource) {
                    (this->setResources_(subResource, nullOrResourcePtr, patch), ...);
                },
                subResources);
        }
    }




    template<typename ResourcesUser, typename ResourcesType>
    void registerResources_(ResourcesUser const& user)
    {
        auto notInMap = [](auto& key, auto& map) {
            auto foundResource = map.find(key);
            return foundResource == std::end(map);
        };

        if constexpr (isUserFieldType<GridLayoutT, ResourcesUser, ResourcesType>::value)
        {
            auto const& resourcesProperties = user.getFieldNamesAndQuantities();
            for (auto const& properties : resourcesProperties)
            {
                std::string const& resourcesName = properties.name;
                auto const& qty                  = properties.qty;

                if (notInMap(resourcesName, nameToResourceInfo_))
                {
                    ResourcesInfo info;
                    info.variable = std::make_shared<typename ResourcesType::variable_type>(
                        resourcesName, qty);

                    info.id = variableDatabase_->registerVariableAndContext(
                        info.variable, context_, SAMRAI::hier::IntVector::getZero(dimension_));

                    nameToResourceInfo_.emplace(resourcesName, info);
                }
            }
        }

        if constexpr (isUserParticleType<ResourcesUser, ResourcesType>::value)
        {
            auto const& resourcesProperties = user.getParticleArrayNames();
            for (auto const& properties : resourcesProperties)
            {
                auto const& name = properties.name;

                if (notInMap(name, nameToResourceInfo_))
                {
                    ResourcesInfo info;

                    info.variable = std::make_shared<typename ResourcesType::variable_type>(name);

                    info.id = variableDatabase_->registerVariableAndContext(
                        info.variable, context_, SAMRAI::hier::IntVector::getZero(dimension_));

                    nameToResourceInfo_.emplace(name, info);
                }
            }
        }
    }



    /** \brief setResourcesInternal_ aims at setting ResourcesUser pointers to the
     * appropriate data on the Patch or to reset them to nullptr.
     */
    template<typename ResourcesUser, typename ResourcesType, typename ResourcesProperties,
             typename RequestedPtr>
    void setResourcesInternal_(ResourcesUser& obj, ResourcesType resourceType,
                               ResourcesProperties const& resourcesProperties,
                               SAMRAI::hier::Patch const& patch, RequestedPtr) const
    {
        for (auto const& properties : resourcesProperties)
        {
            std::string const& resourcesName = properties.name;
            auto const& resourceInfoIt       = nameToResourceInfo_.find(resourcesName);
            if (resourceInfoIt != nameToResourceInfo_.end())
            {
                auto data = getResourcesPointer_<ResourcesType, RequestedPtr>(
                    resourceType, resourceInfoIt->second, patch);

                obj.setBuffer(resourcesName, data);
            }
            else
            {
                throw std::runtime_error("Resources not found !");
            }
        }
    }




    //! \brief Allocate the data on the given level
    template<typename ResourcesUser, typename ResourcesProperties>
    void allocate_(ResourcesUser const& obj, ResourcesProperties const& resourcesProperties,
                   SAMRAI::hier::Patch& patch, double const allocateTime) const
    {
        for (auto const& properties : resourcesProperties)
        {
            std::string const& resourcesName  = properties.name;
            auto const& resourceVariablesInfo = nameToResourceInfo_.find(resourcesName);
            if (resourceVariablesInfo != nameToResourceInfo_.end())
            {
                patch.allocatePatchData(resourceVariablesInfo->second.id, allocateTime);
            }
            else
            {
                throw std::runtime_error("Resources not found !");
            }
        }
    }

    std::string contextName_{"default"};
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    SAMRAI::tbox::Dimension dimension_;
    std::map<std::string, ResourcesInfo> nameToResourceInfo_;

    template<typename ResourcesManager, typename... ResourcesUsers>
    friend class ResourcesGuard;
};

} // namespace PHARE
#endif
