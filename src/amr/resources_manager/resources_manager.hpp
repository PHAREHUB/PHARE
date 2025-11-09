#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_HPP

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include "core/utilities/types.hpp"
#include "field_resource.hpp"
#include "resources_guards.hpp"
#include "particle_resource.hpp"
#include "tensor_field_resource.hpp"
#include "resources_manager_utilities.hpp"

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include "SAMRAI/hier/PatchDataRestartManager.h"

#include <map>
#include <optional>


namespace PHARE
{
namespace amr
{
    /**
     * \brief ResourcesInfo contains the SAMRAI variable and patchdata ID to retrieve data from the
     * SAMRAI database
     */
    struct ResourcesInfo
    {
        std::shared_ptr<SAMRAI::hier::Variable> variable;
        int id;
    };


    struct ResourcesManagerGlobals
    {
        static ResourcesManagerGlobals& INSTANCE();


        std::vector<std::map<std::string, ResourcesInfo>*> resources_;

        static auto ALL_IDS()
        {
            std::vector<int> ids;
            ids.reserve((core::sum_from(INSTANCE().resources_,
                                        [](auto const& map) { return map->size(); })));

            assert(INSTANCE().resources_.size());
            for (auto const res_map_ptr : INSTANCE().resources_)
                for (auto const& [_, info] : *res_map_ptr)
                    ids.emplace_back(info.id);

            return ids;
        }

        static void registerForRestarts()
        {
            auto pdrm = SAMRAI::hier::PatchDataRestartManager::getManager();
            for (auto const& id : ALL_IDS()) // duplicates don't matter
                pdrm->registerPatchDataForRestart(id);
        }
    };


    /** \brief ResourcesManager is an adapter between PHARE objects that manipulate
     * data on patches, and the SAMRAI variable database system, storing the data.
     * It is used by PHARE to register data to the samrai system, allocate data on patches, and to
     * get access to it whenever needed, without having to know the SAMRAI database system.
     *
     * Objects registering and retrieving data through the ResourcesManager are called
     * ResourcesViews. A ResourcesView needs to satisfy a specific interface to work with
     * the ResourcesManager.
     *
     * There are only two kinds of Resources that can be registered to the SAMRAI system via the
     * ResourcesManager, because them only are associated with a PatchData generalization:
     *
     * - Field
     * - ParticleArray
     *
     * Several kinds of ResourcesView can register their resources to the ResourcesManager
     * and are identified by the following type traits:
     *
     * - has_runtime_subresourceview_list <ResourcesView>
     * - has_compiletime_subresourcesview_list <ResourcesView>
     *
     *
     * for example:
     *
     * - has_compiletime_subresourcesview_list<Electromag> is true because it holds 2 VecFields that
     * hold Field objects.
     *
     * ResourcesManager is used to register ResourceUsers to the SAMRAI system, this is done
     * by calling the registerResources() method. It is also used to allocate already registered
     * ResourcesViews on a patch by calling the method allocate(). One can also get the identifier
     * (ID) of a patchdata corresponding to one or several ResourcesViews by calling getID() or
     * getIDs() methods.
     *
     * Data objects like VecField and Ions etc., i.e. ResourcesViews, need to be set on a patch
     * before being usable. Having a Patch and several ResourcesViews obj1, obj2, etc. this is done
     * by calling:
     *
     * dataOnPatch = setOnPatch(patch, obj1, obj2);
     *
     * obj1 and obj2 become unusable again at the end of the scope of dataOnPatch
     *
     *
     */

    template<typename GridLayoutT, typename Grid_t>
    class ResourcesManager
    {
        using This = ResourcesManager<GridLayoutT, Grid_t>;
        using QuantityType =
            typename extract_quantity_type<typename Grid_t::physical_quantity_type>::type;

    public:
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;

        using UserField_t = UserFieldType<Grid_t, GridLayoutT>;

        template<typename ResourcesView>
        using UserParticle_t = UserParticleType<ResourcesView, interp_order>;

        template<std::size_t rank>
        using UserTensorField_t = UserTensorFieldType<rank, Grid_t, GridLayoutT, QuantityType>;


        ResourcesManager()
            : variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
            , context_{variableDatabase_->getContext(contextName_)}
            , dimension_{SAMRAI::tbox::Dimension{dimension}}
        {
            ResourcesManagerGlobals::INSTANCE().resources_.emplace_back(&nameToResourceInfo_);
        }
        ~ResourcesManager()
        {
            for (auto& [key, resourcesInfo] : nameToResourceInfo_)
                variableDatabase_->removeVariable(key);

            auto& vec = ResourcesManagerGlobals::INSTANCE().resources_;
            vec.erase(std::remove(vec.begin(), vec.end(), &nameToResourceInfo_), vec.end());
        }


        ResourcesManager(ResourcesManager const&)                   = delete;
        ResourcesManager(ResourcesManager&&)                        = delete;
        ResourcesManager& operator=(ResourcesManager const& source) = delete;
        ResourcesManager& operator=(ResourcesManager&&)             = delete;


        /** @brief registerResources takes a ResourcesView to register its resources
         * to the SAMRAI data system.
         *
         * the function asks the ResourcesView informations (at least the name) of the
         * buffers to register. This information is used to know which kind of
         * SAMRAI patchdata & variable to create, but is also stored to retrieve
         * this data when the ResourcesView is later going to ask for it (via the BufferGuard)
         *
         * The function registers all FieldData for ResourcesView that have Fields, all ParticleData
         * for ResourcesView that have particle Arrays. In case the ResourcesView has sub-resources
         * we ask for them in a tuple, and recursively call registerResources() for all of the
         * unpacked elements
         */
        template<typename ResourcesView>
        void registerResources(ResourcesView& obj)
        {
            if constexpr (is_resource<ResourcesView>::value)
            {
                registerResource_(obj);
            }
            else
            {
                static_assert(has_sub_resources_v<ResourcesView>);

                if constexpr (has_runtime_subresourceview_list<ResourcesView>::value)
                {
                    for (auto& resourcesUser : obj.getRunTimeResourcesViewList())
                    {
                        this->registerResources(resourcesUser);
                    }
                }

                if constexpr (has_compiletime_subresourcesview_list<ResourcesView>::value)
                {
                    // unpack the tuple subResources and apply for each element registerResources()
                    // (recursively)
                    std::apply(
                        [this](auto&... subResource) {
                            (this->registerResources(subResource), ...);
                        },
                        obj.getCompileTimeResourcesViewList());
                }
            }
        }



        /** @brief allocate the appropriate PatchDatas on the Patch for the ResourcesView
         *
         * The function allocates all FieldData for ResourcesView that have Fields, all ParticleData
         * for ResourcesView that have particle Arrays. In case the ResourcesView has sub-resources
         * we ask for them in a tuple, and recursively call allocate() for all of the unpacked
         * elements
         */
        template<typename ResourcesView>
        void allocate(ResourcesView& obj, SAMRAI::hier::Patch& patch,
                      double const allocateTime) const
        {
            if constexpr (is_resource<ResourcesView>::value)
            {
                allocate_(obj, patch, allocateTime);
            }
            else
            {
                static_assert(has_sub_resources_v<ResourcesView>);

                if constexpr (has_runtime_subresourceview_list<ResourcesView>::value)
                {
                    for (auto& resourcesUser : obj.getRunTimeResourcesViewList())
                    {
                        this->allocate(resourcesUser, patch, allocateTime);
                    }
                }

                if constexpr (has_compiletime_subresourcesview_list<ResourcesView>::value)
                {
                    // unpack the tuple subResources and apply for each element registerResources()
                    std::apply([this, &patch, allocateTime](auto&... subResource) //
                               { (this->allocate(subResource, patch, allocateTime), ...); },
                               obj.getCompileTimeResourcesViewList());
                }
            }
        }




        /** \brief set all passed resources on given Patch
         *
         * auto dataOnPatch = setOnPatch(patch, obj1, obj2, ...);
         *
         * now obj1, obj2 data containers contain data defined on the given patch.
         * At the end of the scope of dataOnPatch, obj1 and obj2 will become unusable again
         */
        template<typename... ResourcesViews>
        NO_DISCARD constexpr ResourcesGuard<ResourcesManager, ResourcesViews...>
        setOnPatch(SAMRAI::hier::Patch const& patch, ResourcesViews&... resourcesUsers)
        {
            return ResourcesGuard<ResourcesManager, ResourcesViews...>{patch, *this,
                                                                       resourcesUsers...};
        }



        /** @brief getTime is used to get the time of the Resources associated with the given
         * ResourcesView on the given patch.
         */
        template<typename ResourcesView>
        NO_DISCARD auto getTimes(ResourcesView& obj, SAMRAI::hier::Patch const& patch) const
        {
            auto IDs = getIDs(obj);
            std::vector<double> times;
            /*std::transform(std::begin(IDs), std::end(IDs), std::back_inserter(times),
                           [&patch](auto const& id) {
                               auto patchdata = patch.getPatchData(id);
                               return patchdata->getTime();
                           });*/

            for (auto const& id : IDs)
            {
                auto patchdata = patch.getPatchData(id);
                times.push_back(patchdata->getTime());
            }
            return times;
        }



        /** @brief setTime is used to set the time of the Resources associated with the given
         * ResourcesView on the given patch.
         */
        template<typename ResourcesView>
        void setTime(ResourcesView& obj, SAMRAI::hier::Patch const& patch, double time) const
        {
            for (auto const& id : getIDs(obj))
            {
                auto patchdata = patch.getPatchData(id);
                patchdata->setTime(time);
            }
        }



        /** \brief Get all the names and resources id that the resource view
         *  have registered via the ResourcesManager
         */
        template<typename ResourcesView>
        NO_DISCARD std::vector<int> getIDs(ResourcesView& obj) const
        {
            std::vector<int> IDs;
            this->getIDs_(obj, IDs);
            return IDs;
        }



        /** \brief Get all the names and resources id that the resource view
         *  have registered via the ResourcesManager
         */
        NO_DISCARD std::optional<int> getID(std::string const& resourceName) const
        {
            auto id = nameToResourceInfo_.find(resourceName);

            if (id != std::end(nameToResourceInfo_))
                return std::optional<int>{id->second.id};

            return std::nullopt;
        }




        void registerForRestarts() const
        {
            auto pdrm = SAMRAI::hier::PatchDataRestartManager::getManager();
            for (auto const& id : restart_patch_data_ids())
                pdrm->registerPatchDataForRestart(id);
        }

        // needed as long as we have different resource managers dealing with different physical
        // quantities
        // template<typename ResourcesView>
        // void registerForRestarts(ResourcesView const& view) const
        // {
        //     auto pdrm = SAMRAI::hier::PatchDataRestartManager::getManager();

        //     for (auto const& id : restart_patch_data_ids(view))
        //         pdrm->registerPatchDataForRestart(id);
        // }


        NO_DISCARD auto restart_patch_data_ids() const
        { // see https://github.com/PHAREHUB/PHARE/issues/664
            return ResourcesManagerGlobals::ALL_IDS();
            // std::vector<int> ids;
            // for (auto const& [key, info] : nameToResourceInfo_)
            //     ids.emplace_back(info.id);
            // return ids;
        }

        // template<typename ResourcesView> // this function is never called
        // NO_DISCARD auto restart_patch_data_ids(ResourcesView const& view) const
        // {
        //     std::vector<int> ids;
        //     getIDs_(view, ids);
        //     return ids;
        // }


        auto getIDsList(auto&&... keys) const
        {
            auto const Fn = [&](auto& key) {
                if (key.empty())
                    throw std::runtime_error("Resource Manager key cannot be empty");
                if (auto const id = getID(key))
                    return *id;
                throw std::runtime_error("Resource Manager has no key: " + key);
            };
            return std::array{Fn(keys)...};
        }


        // iterate per patch and set args on patch
        template<typename... Args>
        auto inline enumerate(SAMRAI::hier::PatchLevel const& level, Args&&... args)
        {
            return LevelLooper<SAMRAI::hier::PatchLevel const, Args...>{*this, level, args...};
        }
        template<typename... Args>
        auto inline enumerate(SAMRAI::hier::PatchLevel& level, Args&&... args)
        {
            return LevelLooper<SAMRAI::hier::PatchLevel, Args...>{*this, level, args...};
        }

    private:
        template<typename Level_t, typename... Args>
        struct LevelLooper
        {
            LevelLooper(ResourcesManager& rm, Level_t& lvl, Args&... arrgs)
                : rm{rm}
                , level{lvl}
                , args{std::forward_as_tuple(arrgs...)}
            {
            }


            struct Iterator
            {
                void operator++() { ++raw; }
                bool operator==(Iterator const& that) { return raw == that.raw; }
                bool operator!=(Iterator const& that) { return raw != that.raw; }
                std::shared_ptr<SAMRAI::hier::Patch> const& operator*()
                {
                    looper->set(**raw);
                    return *raw;
                }

                ~Iterator() { looper->unset(); }

                LevelLooper* looper;
                SAMRAI::hier::PatchLevel::Iterator raw;
            };

            void set(auto& patch)
            {
                std::apply([&](auto&... user) { ((rm.setResources_(user, patch)), ...); }, args);
            }

            void unset()
            {
                std::apply([&](auto&... user) { ((rm.unsetResources_(user)), ...); }, args);
            }

            auto begin() { return Iterator{this, level.begin()}; }
            auto end() { return Iterator{this, level.end()}; };

            ResourcesManager& rm;
            Level_t& level;
            std::tuple<Args&...> args;
        };



        template<typename ResourcesView>
        void getIDs_(ResourcesView& obj, std::vector<int>& IDs) const
        {
            if constexpr (is_resource<ResourcesView>::value)
            {
                auto foundIt = nameToResourceInfo_.find(obj.name());
                if (foundIt == nameToResourceInfo_.end())
                    throw std::runtime_error("Cannot find " + obj.name());
                IDs.push_back(foundIt->second.id);
            }
            else
            {
                if constexpr (has_runtime_subresourceview_list<ResourcesView>::value)
                {
                    for (auto& resourcesUser : obj.getRunTimeResourcesViewList())
                    {
                        //
                        this->getIDs_(resourcesUser, IDs);
                    }
                }

                if constexpr (has_compiletime_subresourcesview_list<ResourcesView>::value)
                {
                    // unpack the tuple subResources and apply for each element registerResources()
                    std::apply(
                        [this, &IDs](auto&... subResource) {
                            (this->getIDs_(subResource, IDs), ...);
                        },
                        obj.getCompileTimeResourcesViewList());
                }
            }
        }


        // The function getResourcesPointer_ is the one that depending
        // on NullOrResourcePtr will choose to return
        // the real pointer or a nullptr to the correct type.

        /** \brief Returns a pointer to the patch data instantiated in the patch
         * is used by getResourcesPointer_ when view code wants the pointer to the data
         */
        template<typename ResourceType>
        auto getPatchData_(ResourcesInfo const& resourcesVariableInfo,
                           SAMRAI::hier::Patch const& patch) const
        {
            auto patchData = patch.getPatchData(resourcesVariableInfo.variable, context_);
            return (std::dynamic_pointer_cast<typename ResourceType::patch_data_type>(patchData))
                ->getPointer();
        }




        /** \brief this specialization of getResourcesPointer_ is used when
         * the client code wants to get a pointer to a patch data resource
         */
        template<typename ResourceType>
        auto getResourcesPointer_(ResourcesInfo const& resourcesVariableInfo,
                                  SAMRAI::hier::Patch const& patch) const
        {
            return getPatchData_<ResourceType>(resourcesVariableInfo, patch);
        }



        void static handle_sub_resources(auto fn, auto& obj, auto&&... args)
        {
            using ResourcesView = decltype(obj);

            if constexpr (has_runtime_subresourceview_list<ResourcesView>::value)
                for (auto& runtimeResource : obj.getRunTimeResourcesViewList())
                    fn(runtimeResource, args...);

            // unpack the tuple subResources and apply for each element registerResources()
            //  (recursively)
            if constexpr (has_compiletime_subresourcesview_list<ResourcesView>::value)
                std::apply([&](auto&... subResource) { (fn(subResource, args...), ...); },
                           obj.getCompileTimeResourcesViewList());
        }


        template<typename ResourcesView>
        void setResources_(ResourcesView& obj, SAMRAI::hier::Patch const& patch) const
        {
            if constexpr (is_resource<ResourcesView>::value)
            {
                setResourcesInternal_(obj, patch);
            }
            else
            {
                static_assert(has_sub_resources_v<ResourcesView>);

                handle_sub_resources( //
                    [&](auto&&... args) { this->setResources_(args...); }, obj, patch);
            }
        }

        template<typename ResourcesView>
        void unsetResources_(ResourcesView& obj) const
        {
            if constexpr (is_resource<ResourcesView>::value)
            {
                unsetResourcesInternal_(obj);
            }
            else
            {
                static_assert(has_sub_resources_v<ResourcesView>);
                handle_sub_resources([&](auto&&... args) { this->unsetResources_(args...); }, obj);
            }
        }


        template<typename ResourcesView>
        void registerResource_(ResourcesView const& view)
        {
            using ResourcesResolver_t = ResourceResolver<This, ResourcesView>;

            if (view.name().empty())
                throw std::runtime_error("Resource Manager key cannot be empty");

            if (nameToResourceInfo_.count(view.name()) == 0)
            {
                ResourcesInfo info;
                info.variable = ResourcesResolver_t::make_shared_variable(view);
                info.id       = variableDatabase_->registerVariableAndContext(
                    info.variable, context_, SAMRAI::hier::IntVector::getZero(dimension_));

                nameToResourceInfo_.emplace(view.name(), info);
            }
        }



        /** \brief setResourcesInternal_ aims at setting ResourcesView pointers to the
         * appropriate data on the Patch or to reset them to nullptr.
         */
        template<typename ResourcesView>
        void setResourcesInternal_(ResourcesView& obj, SAMRAI::hier::Patch const& patch) const
        {
            using ResourceResolver_t = ResourceResolver<This, ResourcesView>;
            using ResourcesType      = ResourceResolver_t::type;

            if (obj.name().empty())
                throw std::runtime_error("Resource Manager key cannot be empty");

            auto const& resourceInfoIt = nameToResourceInfo_.find(obj.name());
            if (resourceInfoIt == nameToResourceInfo_.end())
                throw std::runtime_error("Resources not found ! " + obj.name());

            obj.setBuffer(getResourcesPointer_<ResourcesType>(resourceInfoIt->second, patch));
        }

        template<typename ResourcesView>
        void unsetResourcesInternal_(ResourcesView& obj) const
        {
            if (obj.name().empty())
                throw std::runtime_error("Resource Manager key cannot be empty");

            auto const& resourceInfoIt = nameToResourceInfo_.find(obj.name());
            if (resourceInfoIt == nameToResourceInfo_.end())
                throw std::runtime_error("Resources not found !");

            obj.setBuffer(nullptr);
        }



        //! \brief Allocate the data on the given level
        template<typename ResourcesView>
        void allocate_(ResourcesView const& obj, SAMRAI::hier::Patch& patch,
                       double const allocateTime) const
        {
            std::string const& resourcesName = obj.name();

            if (obj.name().empty())
                throw std::runtime_error("Resource Manager key cannot be empty");

            auto const& resourceVariablesInfo = nameToResourceInfo_.find(resourcesName);
            if (resourceVariablesInfo != nameToResourceInfo_.end())
            {
                if (!patch.checkAllocated(resourceVariablesInfo->second.id))
                    patch.allocatePatchData(resourceVariablesInfo->second.id, allocateTime);
            }
            else
            {
                throw std::runtime_error("Resources not found ! " + resourcesName);
            }
        }

        std::string contextName_{"default"};
        SAMRAI::hier::VariableDatabase* variableDatabase_;
        std::shared_ptr<SAMRAI::hier::VariableContext> context_;
        SAMRAI::tbox::Dimension dimension_;
        std::map<std::string, ResourcesInfo> nameToResourceInfo_;

        template<typename ResourcesManager, typename... ResourcesViews>
        friend class ResourcesGuard;
    };
} // namespace amr
} // namespace PHARE
#endif
