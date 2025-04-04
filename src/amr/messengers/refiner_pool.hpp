#ifndef PHARE_REFINER_POOL_HPP
#define PHARE_REFINER_POOL_HPP

#include "refiner.hpp"


#include <map>
#include <memory>
#include <string>
#include <vector>



namespace PHARE
{
namespace amr
{

    /**
     * A RefinerPool is a container of Refiner objects that can (but not necessarily)
     * be processed together.
     */
    template<typename ResourcesManager, RefinerType Type>
    class RefinerPool
    {
        using Refiner_t      = Refiner<ResourcesManager, Type>;
        using RefineOperator = SAMRAI::hier::RefineOperator;

    public:
        /*@brief add a static communication between sources and destinations.
         * This overload takes several sources/destinations/keys and add one refiner for each*/
        template<typename Names>
        void addStaticRefiners(Names const& destinations, Names const& sources,
                               std::shared_ptr<RefineOperator> refineOp,
                               std::vector<std::string> keys);


        /*@brief convenience overload of the above when source = destination, for VecField*/
        template<typename Names>
        void addStaticRefiners(Names const& src_dest, std::shared_ptr<RefineOperator> refineOp,
                               std::vector<std::string> key);

        /* @brief add a static communication between a single source and destination.*/
        template<typename Name>
        void addStaticRefiner(Name const& ghostName, Name const& src,
                              std::shared_ptr<RefineOperator> const& refineOp,
                              std::string const key);

        /**
         * @brief convenience overload of above addStaticRefiner taking only one name
         * used for communications from a quantity to the same quantity.*/
        template<typename Name>
        void addStaticRefiner(Name const& src_dest, std::shared_ptr<RefineOperator> const& refineOp,
                              std::string const key);


        /**
         * @brief fill the given pool of refiners with a new refiner per VecField
         * in ghostVecs. Data will be spatially refined using the specified refinement
         * operator, and time interpolated between time n and n+1 of next coarser data,
         * represented by modelVec and oldModelVec.*/
        void addTimeRefiners(std::vector<core::VecFieldNames> const& ghostVecs,
                             core::VecFieldNames const& modelVec,
                             core::VecFieldNames const& oldModelVec,
                             std::shared_ptr<RefineOperator>& refineOp,
                             std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>& timeOp);



        /**
         * add a refiner that will use time and spatial interpolation.
         * time interpolation will be done between data represented by model and oldModel
         * , and use the timeOp operator. Spatial refinement of the result
         * will be done using the refineOp operator and the result put in the data
         * represented by `ghost`.
         * The refiner added to the pool will be retrievable using the given key.
         *
         * This overload is for vector fields*/
        void addTimeRefiner(core::VecFieldNames const& ghost, core::VecFieldNames const& model,
                            core::VecFieldNames const& oldModel,
                            std::shared_ptr<RefineOperator> const& refineOp,
                            std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp,
                            std::string key);


        // mhd needs to expose several scalars to the time refiners
        void addTimeRefiners(std::vector<std::string> const& ghostField,
                             std::string const& modelField, std::string const& oldModelField,
                             std::shared_ptr<RefineOperator> const& refineOp,
                             std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp);

        // this overload takes simple strings.
        void addTimeRefiner(std::string const& ghost, std::string const& model,
                            std::string const& oldModel,
                            std::shared_ptr<RefineOperator> const& refineOp,
                            std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp,
                            std::string key);



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
        {
            for (auto& [_, refiner] : refiners_)
            {
                refiner.registerLevel(hierarchy, level);
            }
        }



        /** @brief this overload will execute communications for all quantities in the pool. */
        void fill(int const levelNumber, double const initDataTime) const
        {
            for (auto const& [key, refiner] : refiners_)
            {
                refiner.fill(levelNumber, initDataTime);
            }
        }

        /** @brief filldGhosts is used to fill the ghost nodes of all the given VecField.*/
        template<typename VecFieldT>
        void fill(VecFieldT& vec, int const levelNumber, double const fillTime)
        {
            if (refiners_.count(vec.name()) == 0)
                throw std::runtime_error("no refiner for " + vec.name());

            refiners_.at(vec.name()).fill(vec, levelNumber, fillTime);
        }



        /** @brief executes a regridding for all quantities in the pool.*/
        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            int const levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            double const initDataTime)
        {
            for (auto& [key, refiner] : refiners_)
            {
                refiner.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            }
        }


        RefinerPool(std::shared_ptr<ResourcesManager> const& rm)
            : rm_{rm}
        {
        }


    private:
        using Qty = std::string;
        std::map<Qty, Refiner<ResourcesManager, Type>> refiners_;
        std::shared_ptr<ResourcesManager> rm_;
    };


    template<typename ResourcesManager, RefinerType Type>
    template<typename Names>
    void RefinerPool<ResourcesManager, Type>::addStaticRefiners(
        Names const& destinations, Names const& sources, std::shared_ptr<RefineOperator> refineOp,
        std::vector<std::string> keys)
    {
        assert(destinations.size() == sources.size());
        auto key = std::begin(keys);
        for (std::size_t i = 0; i < destinations.size(); ++i)
        {
            addStaticRefiner(destinations[i], sources[i], refineOp, *key++);
        }
    }


    template<typename ResourcesManager, RefinerType Type>
    template<typename Names>
    void
    RefinerPool<ResourcesManager, Type>::addStaticRefiners(Names const& src_dest,
                                                           std::shared_ptr<RefineOperator> refineOp,
                                                           std::vector<std::string> key)
    {
        addStaticRefiners(src_dest, src_dest, refineOp, key);
    }




    template<typename ResourcesManager, RefinerType Type>
    template<typename Name>
    void RefinerPool<ResourcesManager, Type>::addStaticRefiner(
        Name const& ghostName, Name const& src, std::shared_ptr<RefineOperator> const& refineOp,
        std::string const key)
    {
        auto const [it, success]
            = refiners_.insert({key, Refiner_t(ghostName, src, rm_, refineOp)});

        if (!success)
            throw std::runtime_error(key + " is already registered");
    }



    template<typename ResourcesManager, RefinerType Type>
    template<typename Name>
    void RefinerPool<ResourcesManager, Type>::addStaticRefiner(
        Name const& descriptor, std::shared_ptr<RefineOperator> const& refineOp,
        std::string const key)
    {
        addStaticRefiner(descriptor, descriptor, refineOp, key);
    }


    template<typename ResourcesManager, RefinerType Type>
    void RefinerPool<ResourcesManager, Type>::addTimeRefiners(
        std::vector<core::VecFieldNames> const& ghostVecs, core::VecFieldNames const& modelVec,
        core::VecFieldNames const& oldModelVec, std::shared_ptr<RefineOperator>& refineOp,
        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>& timeOp)
    {
        for (auto const& ghostVec : ghostVecs)
        {
            addTimeRefiner(ghostVec, modelVec, oldModelVec, refineOp, timeOp, ghostVec.vecName);
        }
    }

    template<typename ResourcesManager, RefinerType Type>
    void RefinerPool<ResourcesManager, Type>::addTimeRefiner(
        core::VecFieldNames const& ghost, core::VecFieldNames const& model,
        core::VecFieldNames const& oldModel, std::shared_ptr<RefineOperator> const& refineOp,
        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp, std::string key)
    {
        auto const [it, success]
            = refiners_.insert({key, Refiner_t(ghost, model, oldModel, rm_, refineOp, timeOp)});
        if (!success)
            throw std::runtime_error(key + " is already registered");
    }

    template<typename ResourcesManager, RefinerType Type>
    void RefinerPool<ResourcesManager, Type>::addTimeRefiners(
        std::vector<std::string> const& ghostFields, std::string const& modelField,
        std::string const& oldModelField, std::shared_ptr<RefineOperator> const& refineOp,
        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp)
    {
        for (auto const& ghostField : ghostFields)
        {
            addTimeRefiner(ghostField, modelField, oldModelField, refineOp, timeOp, ghostField);
        }
    }

    template<typename ResourcesManager, RefinerType Type>
    void RefinerPool<ResourcesManager, Type>::addTimeRefiner(
        std::string const& ghost, std::string const& model, std::string const& oldModel,
        std::shared_ptr<RefineOperator> const& refineOp,
        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp, std::string key)
    {
        auto const [it, success]
            = refiners_.insert({key, Refiner_t(ghost, model, oldModel, rm_, refineOp, timeOp)});
        if (!success)
            throw std::runtime_error(key + " is already registered");
    }
} // namespace amr
} // namespace PHARE

#endif
