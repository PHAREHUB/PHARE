#ifndef PHARE_AMR_TOOLS_RESOURCES_MANAGER_HPP
#define PHARE_AMR_TOOLS_RESOURCES_MANAGER_HPP


#include "field_resource.hpp"
#include "particle_resource.hpp"
#include "resources_guards.hpp"
#include "resources_manager_utilities.hpp"
#include "core/def.hpp"


#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/VariableDatabase.h>
#include "SAMRAI/hier/PatchDataRestartManager.h"

#include <map>
#include <optional>
#include <variant>


namespace PHARE
{
namespace amr
{

    template<typename FieldView, typename... OtherViews>
    class HierarchyViews
    {
        using AnyView = std::variant<FieldView, OtherViews...>;

    public:
        using PatchViews = std::vector<AnyView>;
        using LevelViews = std::vector<PatchViews>;
        void add(std::string const& name, std::size_t const levelNumber, auto viewPtr)
        {
            assert(exists_(name));
            if (isNewLevel_(views_[name], levelNumber))
            {
                views_[name].emplace_back(PatchViews{});
            }
            views_[name][levelNumber].emplace_back(*viewPtr);
        }

        // \brief deletes all views for a given level
        void resetLevel(std::string const& name, std::size_t const levelNumber)
        {
            assert(exists_(name));
            if (levelExists_(views_[name], levelNumber))
                views_[name][levelNumber].clear();
        }

    private:
        bool exists_(std::string const& name) const
        {
            return views_.find(name) != std::end(views_);
        }

        bool isNewLevel_(LevelViews const& lvlViews, std::size_t const ilvl) const
        {
            assert(ilvl <= lvlViews.size());
            return ilvl == views_.size();
        }

        bool levelExists_(LevelViews const& lvlViews, std::size_t const ilvl) const
        {
            return !isNewLevel_(lvlViews, ilvl);
        }

        std::unordered_map<std::string, LevelViews> views_;
    };


    /** \brief ResourcesManager is an adapter between PHARE objects that manipulate
     * data on patches, and the SAMRAI variable database system, storing the data.
     * It is used by PHARE to register data to the samrai system, allocate data on patches, and
     * to get access to it whenever needed, without having to know the SAMRAI database system.
     *
     * Objects registering and retrieving data through the ResourcesManager are called
     * Views. Views give PHARE access to data stored on patches, called Resources.
     *
     * There are only two kinds of Resources that can be registered to the SAMRAI system via the
     * ResourcesManager, because them only are associated with a PatchData generalization:
     *
     * - Grid, allocates meshed data for the FieldData
     * - ParticlePack, holds particle allocated in ParticleArrays in the ParticleData
     *
     * A View needs to satisfy a specific interface to work with the ResourcesManager.
     * Several kinds of View can register their resources to the ResourcesManager
     * and are identified by the following type traits:
     *
     * - has_runtime_subresourceview_list <ResourcesView>
     * - has_compiletime_subresourcesview_list <ResourcesView>
     *
     *
     * for example:
     *
     * - has_compiletime_subresourcesview_list<Electromag> is true because it holds 2 VecFields
     * that hold Field objects.
     *
     * ResourcesManager is used to register Resources to the SAMRAI system, this is done
     * by passing Views attached to these resources to the registerResources() method.
     * It is also used to allocate already registered Resources on a patch by calling the method
     * allocate(). One can also get the identifier (ID) of a patchdata corresponding to one or
     * several Resources by calling getID() or getIDs() methods.
     *
     * Data objects like VecField and Ions etc., i.e. ResourcesViews, need to be set on a patch
     * before being usable. Having a Patch and several ResourcesViews obj1, obj2, etc. this is
     * done by calling:
     *
     * dataOnPatch = setOnPatch(patch, obj1, obj2);
     *
     * obj1 and obj2 become unusable again at the end of the scope of dataOnPatch
     *
     *
     */

    class IResourcesManager
    {
    public:
        virtual void makeLevelViews(SAMRAI::hier::PatchHierarchy const& hierarchy,
                                    int const levelNumber)
            = 0;

        virtual ~IResourcesManager() = default;
    };


    template<typename GridLayoutT, typename Grid_t, typename ParticleArray>
    struct ViewMetaWithParticles
    {
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;

        using FieldViewInfo_t    = FieldViewInfo<Grid_t, GridLayoutT>;
        using ParticleViewInfo_t = ParticleViewInfo<ParticleArray, interp_order>;

        using FieldView_t                   = typename FieldViewInfo_t::view_type;
        using ParticleView_t                = typename ParticleViewInfo_t::view_type;
        static bool constexpr has_particles = true;

        using HierarchyViews_t = HierarchyViews<FieldView_t, ParticleView_t>;
    };


    template<typename GridLayoutT, typename Grid_t>
    struct ViewMetaMeshOnly
    {
        static constexpr std::size_t dimension    = GridLayoutT::dimension;
        static constexpr std::size_t interp_order = GridLayoutT::interp_order;

        using FieldViewInfo_t               = FieldViewInfo<Grid_t, GridLayoutT>;
        using FieldView_t                   = typename FieldViewInfo_t::view_type;
        static bool constexpr has_particles = false;

        using HierarchyViews_t = HierarchyViews<FieldView_t>;
    };


    template<typename GridLayoutT, typename Grid_t, typename ParticleArray = void>
    using SuperDispatcher
        = std::conditional_t<std::is_same_v<ParticleArray, void>,
                             ViewMetaMeshOnly<GridLayoutT, Grid_t>,
                             ViewMetaWithParticles<GridLayoutT, Grid_t, ParticleArray>>;

    template<typename GridLayoutT, typename Grid_t, typename ParticleArray = void>
    class ResourcesManager : public SuperDispatcher<GridLayoutT, Grid_t, ParticleArray>,
                             public IResourcesManager
    {
        using This = ResourcesManager<GridLayoutT, Grid_t, ParticleArray>;

    public:
        using Super                     = SuperDispatcher<GridLayoutT, Grid_t, ParticleArray>;
        static constexpr auto dimension = Super::dimension;
        static constexpr auto interp    = Super::interp_order;

        enum class Resource_t { Field, Particle };

        /**
         * \brief ResourcesInfo contains the SAMRAI variable and patchdata ID to retrieve data
         * from the SAMRAI database
         */
        struct ResourcesInfo
        {
            std::shared_ptr<SAMRAI::hier::Variable> variable;
            int id;
            Resource_t type;
        };

        ResourcesManager()
            : variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
            , context_{variableDatabase_->getContext(contextName_)}
            , dimension_{SAMRAI::tbox::Dimension{dimension}}
        {
        }


        ResourcesManager(ResourcesManager const&)                   = delete;
        ResourcesManager(ResourcesManager&&)                        = delete;
        ResourcesManager& operator=(ResourcesManager const& source) = delete;
        ResourcesManager& operator=(ResourcesManager&&)             = delete;

        template<typename Fn, typename... Args>
        auto recursive_on_view(auto& view, Fn& fn, Args&... args)
        {
            if constexpr (is_final_view_v<decltype(view)>)
            {
                fn(view, std::forward<Args>(args)...);
            }
            else
            {
                if constexpr (has_runtime_views_v<decltype(view)>)
                {
                    for (auto& subview : view.getRunTimeResourcesViewList())
                    {
                        recursive_on_view(subview, fn, args...);
                    }
                }
                if constexpr (has_compiletime_views_v<decltype(view)>)
                {
                    std::apply(
                        [&](auto&... subviews) { (recursive_on_view(subviews, fn, args...), ...); },
                        view.getCompileTimeResourcesViewList());
                }
            }
        }

        /** @brief registerResources takes a View to register its resources
         * to the SAMRAI data system.
         *
         * the function asks the View informations (at least the name) of the
         * resources to register. This information is used to know which kind of
         * SAMRAI patchdata & variable to create, but is also stored to retrieve
         * this data when the View is later going to ask for it.
         */
        template<typename View>
        void registerResources(View& view)
        {
            if constexpr (is_final_view<View>::value)
            {
                registerResource_(view);
            }
            else
            {
                static_assert(has_subview_v<View>);

                if constexpr (has_runtime_views_v<View>)
                {
                    for (auto& resourcesUser : view.getRunTimeResourcesViewList())
                    {
                        this->registerResources(resourcesUser);
                    }
                }

                if constexpr (has_compiletime_views_v<View>)
                {
                    // unpack the tuple subResources and apply for each element
                    // registerResources() (recursively)
                    std::apply([this](auto&... views) { (this->registerResources(views), ...); },
                               view.getCompileTimeResourcesViewList());
                }
            }
        }



        /** @brief allocate the appropriate PatchDatas on the Patch for the View
         *
         * The function allocates all FieldData for View that have Fields, all ParticleData
         * for a View that have particle packs. In case the View has subviews
         * we ask for them in a tuple, and recursively call allocate()          */
        template<typename View>
        void allocate(View& view, SAMRAI::hier::Patch& patch, double const allocateTime) const
        {
            if constexpr (is_final_view_v<View>)
            {
                allocate_(view, patch, allocateTime);
            }
            else
            {
                static_assert(has_subview_v<View>);

                if constexpr (has_runtime_views_v<View>)
                {
                    for (auto& subview : view.getRunTimeResourcesViewList())
                    {
                        this->allocate(subview, patch, allocateTime);
                    }
                }

                if constexpr (has_compiletime_views_v<View>)
                {
                    // unpack the tuple subResources and apply for each element
                    // registerResources()
                    std::apply([this, &patch, allocateTime](auto&... subviews) //
                               { (this->allocate(subviews, patch, allocateTime), ...); },
                               view.getCompileTimeResourcesViewList());
                }
            }
        }




        /** \brief attach passed view on resources on the given Patch
         *
         * auto onPatchGuard = setOnPatch(patch, view1, view2, ...);
         *
         * now view1 and view2 and all their possible subviews are attached to the resources of
         * that patch At the end of the scope of onPatchGuard, view1 and view2 will become
         * unusable again
         */
        template<typename... Views>
        NO_DISCARD constexpr ViewGuard<ResourcesManager, Views...>
        setOnPatch(SAMRAI::hier::Patch const& patch, Views&... views)
        {
            return ViewGuard<ResourcesManager, Views...>{patch, *this, views...};
        }



        /** @brief getTime is used to get the time of the Resources associated with the given
         * View on the given patch.
         */
        template<typename View>
        NO_DISCARD auto getTimes(View& obj, SAMRAI::hier::Patch const& patch) const
        {
            auto IDs = getIDs(obj);
            std::vector<double> times;

            for (auto const& id : IDs)
            {
                auto patchdata = patch.getPatchData(id);
                times.push_back(patchdata->getTime());
            }
            return times;
        }



        /** @brief setTime is used to set the time of the Resources associated with the given
         * View on the given patch.
         */
        template<typename View>
        void setTime(View& obj, SAMRAI::hier::Patch const& patch, double time) const
        {
            for (auto const& id : getIDs(obj))
            {
                auto patchdata = patch.getPatchData(id);
                patchdata->setTime(time);
            }
        }


        template<typename... Views>
        void setTime(SAMRAI::hier::PatchLevel const& level, double time, Views&... views)
        {
            // we take a variadic number of Views, each may represent a runtime/compiletime
            // hierachy of views.
            auto time_setter = [&](auto& view, auto& level, auto time) {
                for (auto& patch : *level)
                    this->setTime(view, *patch, time);
            };
            (recursive_on_view(views, time_setter, level, time), ...);
        }


        /** \brief Get all the names and resources id that the view
         *  have registered via the ResourcesManager
         */
        template<typename View>
        NO_DISCARD std::vector<int> getIDs(View& obj) const
        {
            std::vector<int> IDs;
            this->getIDs_(obj, IDs);
            return IDs;
        }



        /** \brief Get all the names and resources id that the view
         *  have registered via the ResourcesManager
         */
        NO_DISCARD std::optional<int> getID(std::string const& viewName) const
        {
            auto id = nameToResourceInfo_.find(viewName);

            if (id != std::end(nameToResourceInfo_))
                return std::optional<int>{id->second.id};

            return std::nullopt;
        }



        ~ResourcesManager()
        {
            for (auto& [key, resourcesInfo] : nameToResourceInfo_)
            {
                variableDatabase_->removeVariable(key);
            }
        }


        template<typename View>
        void registerForRestarts(View const& view) const
        {
            auto pdrm = SAMRAI::hier::PatchDataRestartManager::getManager();

            for (auto const& id : restart_patch_data_ids(view))
                pdrm->registerPatchDataForRestart(id);
        }

        template<typename View>
        NO_DISCARD auto restart_patch_data_ids(View const& view) const
        {
            // true for now with https://github.com/PHAREHUB/PHARE/issues/664
            constexpr bool ALL_IDS = true;

            std::vector<int> ids;

            if constexpr (ALL_IDS)
            { // get all registered ids to save
                for (auto const& [key, info] : nameToResourceInfo_)
                    ids.emplace_back(info.id);
            }
            else
            { // this is the case when transient datas not to be saved
                getIDs_(view, ids);
            }
            return ids;
        }



        void makeAllViews(SAMRAI::hier::PatchHierarchy const& hierarchy)
        {
            auto const numberOfLevels = hierarchy.getNumberOfLevels();
            for (auto const& [res_name, res_info] : nameToResourceInfo_)
            {
                // TODO Note : this can't work if resoruce manager does not span the
                // whole hierarchy...
                for (auto levelNumber = 0; levelNumber < numberOfLevels; ++levelNumber)
                {
                    makeResourceLevelViews(hierarchy, res_name, res_info, levelNumber);
                }
            }
        }

        void makeLevelViews(SAMRAI::hier::PatchHierarchy const& hierarchy,
                            int const levelNumber) override
        {
            for (auto const& [res_name, res_info] : nameToResourceInfo_)
            {
                makeResourceLevelViews(hierarchy, res_name, res_info, levelNumber);
            }
        }


        void makeResourceLevelViews(SAMRAI::hier::PatchHierarchy const& hierarchy,
                                    std::string const& name, ResourcesInfo const& res_info,
                                    int const levelNumber)
        {
            views_.resetLevel(name, levelNumber);
            auto const patchLevel = hierarchy.getPatchLevel(levelNumber);
            for (auto patch : *patchLevel)
            {
                auto const patchdata_id = res_info.id;
                auto pdata              = patch->getPatchData(patchdata_id);
                // take the name and check if it is a field or a particle
                // if it is a field, we can create a FieldView else a particlepack
                // and set it on the patch
                auto resourceType = res_info.type;
                if (resourceType == Resource_t::Field)
                {
                    auto const viewPtr
                        = getResourcesPtr_<typename Super::FieldViewInfo_t>(res_info, *patch);
                    views_.add(name, levelNumber, viewPtr);
                }
                if constexpr (Super::has_particles)
                    if (resourceType == Resource_t::Particle)
                    {
                        auto const viewPtr = getResourcesPtr_<typename Super::ParticleViewInfo_t>(
                            res_info, *patch);
                        views_.add(name, levelNumber, viewPtr);
                    }
            }
        }


        // return a tuple of references to the levelviews of the given views
        template<typename... Views>
        auto views(SAMRAI::hier::PatchLevel& level, Views&... views)
        {
            auto get_levelview_ref
                = [&](auto const& v) -> typename Super::HierarchyViews_t::PatchViews& {
                auto const& name = v.name();
                return views_.views_[name][level.getLevelNumber()];
            };

            return std::forward_as_tuple(recursive_on_view(views, get_levelview_ref)...);
        }


    private:
        template<typename View>
        void getIDs_(View& view, std::vector<int>& IDs) const
        {
            if constexpr (is_final_view<View>::value)
            {
                auto foundIt = nameToResourceInfo_.find(view.name());
                if (foundIt == nameToResourceInfo_.end())
                    throw std::runtime_error("Cannot find " + view.name());
                IDs.push_back(foundIt->second.id);
            }
            else
            {
                static_assert(has_subview_v<View>);

                if constexpr (has_runtime_views_v<View>)
                {
                    for (auto& views : view.getRunTimeResourcesViewList())
                    {
                        this->getIDs_(views, IDs);
                    }
                }

                if constexpr (has_compiletime_views_v<View>)
                {
                    // unpack the tuple subResources and apply for each element
                    // registerResources()
                    std::apply([this, &IDs](auto&... views) { (this->getIDs_(views, IDs), ...); },
                               view.getCompileTimeResourcesViewList());
                }
            }
        }


        // The function getResourcesPointer_ is the one that depending
        // on NullOrResourcePtr will choose to return
        // the real pointer or a nullptr to the correct type.

        /** \brief Returns a pointer to the patch data instantiated in the patch
         * is used by getResourcesPointer_ when view code wants the pointer to the data
         */
        template<typename ViewInfo>
        auto getResourcesPtr_(ResourcesInfo const& resourcesVariableInfo,
                              SAMRAI::hier::Patch const& patch) const
        {
            auto patchData = patch.getPatchData(resourcesVariableInfo.variable, context_);
            return (std::dynamic_pointer_cast<typename ViewInfo::patch_data_type>(patchData))
                ->getPointer();
        }




        /** \brief returns either the resource pointer or nullptr depending on whether
         * we want to attach or detach a view from the resource
         */
        template<typename ViewInfo, typename RequestedPtr>
        auto getResourceOrNullPtr_(ResourcesInfo const& resourcesVariableInfo,
                                   SAMRAI::hier::Patch const& patch) const
        {
            auto ptr = getResourcesPtr_<ViewInfo>(resourcesVariableInfo, patch);
            auto constexpr attachToResource   = std::is_same_v<RequestedPtr, UseResourcePtr>;
            auto constexpr detachFromResource = std::is_same_v<RequestedPtr, UseNullPtr>;

            if constexpr (attachToResource)
                return ptr;

            else if constexpr (detachFromResource)
                return static_cast<decltype(ptr)>(nullptr);
        }




        template<typename View, typename NullOrResourcePtr>
        void setResources_(View& view, NullOrResourcePtr nullOrResourcePtr,
                           SAMRAI::hier::Patch const& patch) const
        {
            if constexpr (is_final_view_v<View>)
            {
                setResourcesInternal_(view, patch, nullOrResourcePtr);
            }
            else
            {
                static_assert(has_subview_v<View>);

                if constexpr (has_runtime_views_v<View>)
                {
                    for (auto& subview : view.getRunTimeResourcesViewList())
                    {
                        this->setResources_(subview, nullOrResourcePtr, patch);
                    }
                }

                if constexpr (has_compiletime_views_v<View>)
                {
                    std::apply(
                        [this, &patch, &nullOrResourcePtr](auto&... subview) {
                            (this->setResources_(subview, nullOrResourcePtr, patch), ...);
                        },
                        view.getCompileTimeResourcesViewList());
                }
            }
        }


        template<typename View>
        void registerResource_(View const& view)
        {
            using ViewInfoResolver_t = ViewInfoResolver<This, View>;

            if (nameToResourceInfo_.count(view.name()) == 0)
            {
                ResourcesInfo info;
                info.variable = ViewInfoResolver_t::make_shared_variable(view);
                info.id       = variableDatabase_->registerVariableAndContext(
                    info.variable, context_, SAMRAI::hier::IntVector::getZero(dimension_));

                if constexpr (is_field_v<View>)
                {
                    info.type = Resource_t::Field;
                }
                else if constexpr (is_particles_v<View>)
                {
                    info.type = Resource_t::Particle;
                }
                nameToResourceInfo_.emplace(view.name(), info);
            }
        }



        /** \brief setResourcesInternal_ aims at setting ResourcesView pointers to the
         * appropriate data on the Patch or to reset them to nullptr.
         */
        template<typename View, typename RequestedPtr>
        void setResourcesInternal_(View& view, SAMRAI::hier::Patch const& patch, RequestedPtr) const
        {
            using ViewInfo = ViewInfoResolver<This, View>::view_info_type;

            auto const& resourceInfoIt = nameToResourceInfo_.find(view.name());
            if (resourceInfoIt == nameToResourceInfo_.end())
                throw std::runtime_error("Resources not found !");

            view.setBuffer(
                getResourceOrNullPtr_<ViewInfo, RequestedPtr>(resourceInfoIt->second, patch));
        }



        //! \brief Allocate the data on the given level
        template<typename View>
        void allocate_(View const& view, SAMRAI::hier::Patch& patch,
                       double const allocateTime) const
        {
            std::string const& resourcesName  = view.name();
            auto const& resourceVariablesInfo = nameToResourceInfo_.find(resourcesName);
            if (resourceVariablesInfo != nameToResourceInfo_.end())
            {
                if (!patch.checkAllocated(resourceVariablesInfo->second.id))
                    patch.allocatePatchData(resourceVariablesInfo->second.id, allocateTime);
            }
            else
            {
                throw std::runtime_error("Resources not found !");
            }
        }

        std::string contextName_{"default"};
        SAMRAI::hier::VariableDatabase* variableDatabase_;
        std::shared_ptr<SAMRAI::hier::VariableContext> context_;
        SAMRAI::tbox::Dimension dimension_;
        std::map<std::string, ResourcesInfo> nameToResourceInfo_;
        typename Super::HierarchyViews_t views_;

        template<typename ResourcesManager, typename... Views>
        friend class ViewGuard;
    };
} // namespace amr
} // namespace PHARE
#endif
