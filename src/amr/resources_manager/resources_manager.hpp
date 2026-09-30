#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_HPP


#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/def.hpp"
#include "core/utilities/variants.hpp"

#include "amr/samrai.hpp"

#include "resources_guards.hpp"
#include "resources_manager_utilities.hpp"

#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include "SAMRAI/hier/PatchDataRestartManager.h"

#include <map>
#include <array>
#include <tuple>
#include <variant>
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
     * (ID) of a patchdata corresponding to one or several ResourcesViews by calling getIDs(), or
     * by name within one model via scoped<Model>().getID() / getIDsList().
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

    template<auto opts, typename... core_types_per_model>
    class ResourcesManager
    {
        using This = ResourcesManager<opts, core_types_per_model...>;

    public:
        using TypeTuple = std::tuple<core_types_per_model...>;

        /** \brief the UserField_t/UserTensorField_t/UserParticle_t for a single model (e.g.
         * HybridModel/MHDModel) registered on this ResourcesManager - resolved by looking Model
         * up among core_types_per_model, not by re-deriving its shape independently.
         */
        template<typename Model>
        using ResourceUserTypesFor = typename ResourceUserTypesFor_<Model, TypeTuple>::type;

        static constexpr std::size_t dimension = opts.dimension;


        ResourcesManager()
            : variableDatabase_{SamraiLifeCycle::getDatabase()}
            , context_{variableDatabase_->getContext(contextName_)}
            , dimension_{SAMRAI::tbox::Dimension{dimension}}
        {
        }


        ~ResourcesManager()
        {
            for (auto& infos : nameToResourceInfo_)
                for (auto& [_, resourcesInfo] : infos)
                    variableDatabase_->removeVariable(resourcesInfo.variable->getName());
        }


        /** \brief IDs of every resource registered on this ResourcesManager.
         */
        NO_DISCARD auto ALL_IDS() const
        {
            std::vector<int> ids;
            for (auto const& infos : nameToResourceInfo_)
                for (auto const& [_, info] : infos)
                    ids.emplace_back(info.id);

            return ids;
        }

        /** \brief registers, for SAMRAI restarts, the resources of this ResourcesManager.
         */
        void registerForRestarts() const
        {
            auto pdrm = SAMRAI::hier::PatchDataRestartManager::getManager();
            for (auto const& id : ALL_IDS()) // duplicates don't matter
                pdrm->registerPatchDataForRestart(id);
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

        void static handle_sub_resources(auto fn, auto& obj, auto&&... args)
        {
            using ResourcesView = decltype(obj);

            if constexpr (has_runtime_subresourceview_list<ResourcesView>::value)
            {
                for (auto& runtimeResource : obj.getRunTimeResourcesViewList())
                {
                    if constexpr (core::std_variant<decltype(runtimeResource)>)
                        std::visit([&](auto&& val) { fn(val, args...); }, runtimeResource);
                    else
                        fn(runtimeResource, args...);
                }
            }

            if constexpr (has_compiletime_subresourcesview_list<ResourcesView>::value)
            {
                // unpack the tuple subResources and apply for each element registerResources()
                // (recursively)
                std::apply([&](auto&... subResource) { (fn(subResource, args...), ...); },
                           obj.getCompileTimeResourcesViewList());
            }
        }

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

                handle_sub_resources( //
                    [&](auto&&... args) { this->registerResources(args...); }, obj);
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

                handle_sub_resources( //
                    [&](auto&&... args) { this->allocate(args...); }, obj, patch, allocateTime);
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
            std::vector<double> times;
            for (auto const& id : getIDs(obj))
                times.push_back(patch.getPatchData(id)->getTime());
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



        NO_DISCARD auto restart_patch_data_ids() const
        {
            // see https://github.com/PHAREHUB/PHARE/issues/664
            return ALL_IDS();
        }


        // by name within one model, the same name may exist on other models
        template<typename Model>
        NO_DISCARD std::optional<int> getIDFor(std::string const& resourceName) const
        {
            static_assert(model_index_of_<Model> < std::tuple_size_v<TypeTuple>,
                          "Model is not one of the models registered on this ResourcesManager");
            auto const& infos = nameToResourceInfo_[model_index_of_<Model>];
            if (auto it = infos.find(resourceName); it != infos.end())
                return it->second.id;
            return std::nullopt;
        }

        template<typename Model>
        auto getIDsListFor(auto&&... keys) const
        {
            auto const Fn = [&](auto const& key) {
                if (key.empty())
                    throw std::runtime_error("Resource Manager key cannot be empty");
                if (auto const id = getIDFor<Model>(key))
                    return *id;
                throw std::runtime_error("Resource Manager has no key: " + key);
            };
            return std::array{Fn(keys)...};
        }


        /** \brief this ResourcesManager, with name lookups restricted to Model's resources - for
         * callers that only deal with one model (messengers, model views), so names never need
         * to differ between models
         */
        template<typename Model>
        struct ModelScoped
        {
            ResourcesManager const& rm;

            NO_DISCARD auto getID(std::string const& name) const
            {
                return rm.template getIDFor<Model>(name);
            }

            NO_DISCARD auto getIDsList(auto&&... keys) const
            {
                return rm.template getIDsListFor<Model>(keys...);
            }
        };

        template<typename Model>
        NO_DISCARD ModelScoped<Model> scoped() const
        {
            return {*this};
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
            using value_type = std::shared_ptr<SAMRAI::hier::Patch>;

            LevelLooper(ResourcesManager& rm, Level_t& lvl, Args&... arrgs)
                : rm{rm}
                , level{lvl}
                , args{std::forward_as_tuple(arrgs...)}
            {
            }


            struct Patch
            {
                Patch(LevelLooper* looper, std::shared_ptr<SAMRAI::hier::Patch> const& patch)
                    : looper{looper}
                    , patch{patch}
                {
                    looper->set(*patch);
                }

                Patch(Patch const&)            = delete;
                Patch(Patch&&)                 = delete;
                Patch& operator=(Patch const&) = delete;
                Patch& operator=(Patch&&)      = delete;

                SAMRAI::hier::Patch const& operator*() { return *patch; }

                ~Patch() { looper->unset(); }

                LevelLooper* looper;
                std::shared_ptr<SAMRAI::hier::Patch> const patch;
            };

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
            auto size() const { return static_cast<std::size_t>(level.getLocalNumberOfPatches()); }
            auto operator[](std::size_t const idx) { return Patch{this, level.getPatch(idx)}; }

            ResourcesManager& rm;
            Level_t& level;
            std::tuple<Args&...> args;
        };



        template<typename ResourcesView>
        void getIDs_(ResourcesView& obj, std::vector<int>& IDs) const
        {
            if constexpr (is_resource<ResourcesView>::value)
            {
                auto const& infos = resourceInfosFor_(obj);
                auto foundIt      = infos.find(obj.name());
                if (foundIt == infos.end())
                    throw std::runtime_error("Cannot find " + obj.name());
                IDs.push_back(foundIt->second.id);
            }
            else
            {
                static_assert(has_sub_resources_v<ResourcesView>);

                handle_sub_resources([&](auto&&... args) { this->getIDs_(args...); }, obj, IDs);
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

            if (auto& infos = resourceInfosFor_(view); infos.count(view.name()) == 0)
            {
                ResourcesInfo info;
                info.variable = ResourcesResolver_t::make_shared_variable(view);
                info.id       = variableDatabase_->registerVariableAndContext(
                    info.variable, context_, SAMRAI::hier::IntVector::getZero(dimension_));

                infos.emplace(view.name(), info);
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

            auto const& infos          = resourceInfosFor_(obj);
            auto const& resourceInfoIt = infos.find(obj.name());
            if (resourceInfoIt == infos.end())
                throw std::runtime_error("Resources not found ! " + obj.name());

            obj.setBuffer(getResourcesPointer_<ResourcesType>(resourceInfoIt->second, patch));
        }

        template<typename ResourcesView>
        void unsetResourcesInternal_(ResourcesView& obj) const
        {
            if (obj.name().empty())
                throw std::runtime_error("Resource Manager key cannot be empty");

            auto const& infos = resourceInfosFor_(obj);
            if (infos.find(obj.name()) == infos.end())
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

            auto const& infos                 = resourceInfosFor_(obj);
            auto const& resourceVariablesInfo = infos.find(resourcesName);
            if (resourceVariablesInfo != infos.end())
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
        template<typename Model>
        static constexpr std::size_t model_index_of_
            = index_in_tuple<typename model_core_types_of<Model, TypeTuple>::type, TypeTuple>::value;

        // one map per model: names are unique within a model, not across models
        std::array<std::map<std::string, ResourcesInfo>, std::tuple_size_v<TypeTuple>>
            nameToResourceInfo_;

        // the map of the model owning ResourcesView, resolved at compile time
        template<typename ResourcesView>
        NO_DISCARD std::map<std::string, ResourcesInfo>& resourceInfosFor_(ResourcesView const&)
        {
            return nameToResourceInfo_[ResourceResolver<This, ResourcesView>::model_index];
        }

        template<typename ResourcesView>
        NO_DISCARD std::map<std::string, ResourcesInfo> const&
        resourceInfosFor_(ResourcesView const&) const
        {
            return nameToResourceInfo_[ResourceResolver<This, ResourcesView>::model_index];
        }

        template<typename ResourcesManager, typename... ResourcesViews>
        friend class ResourcesGuard;
    };
} // namespace amr
} // namespace PHARE
#endif
