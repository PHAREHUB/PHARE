#ifndef PHARE_REFINER_POOL_HPP
#define PHARE_REFINER_POOL_HPP

#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "refiner.hpp"


#include <map>
#include <memory>
#include <string>
#include <vector>

#include <SAMRAI/xfer/VariableFillPattern.h>


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
        using Refiner_t = Refiner<ResourcesManager, Type>;


    public:
        RefinerPool(std::shared_ptr<ResourcesManager> const& rm)
            : rm_{rm}
        {
        }

        virtual ~RefinerPool() {}
        RefinerPool(RefinerPool const&) = delete;
        RefinerPool(RefinerPool&&)      = default;

        /* @brief add a static communication between a single source and destination.*/
        template<typename Resource, typename Key>
        void
        addStaticRefiner(Resource const& ghostName, Resource const& src,
                         std::shared_ptr<SAMRAI::hier::RefineOperator> const& refineOp,
                         Key const& key,
                         std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
                         std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr);

        /**
         * @brief convenience overload of above addStaticRefiner taking only one name
         * used for communications from a quantity to the same quantity.*/
        template<typename Resource, typename Key>
        void
        addStaticRefiner(Resource const& src_dest,
                         std::shared_ptr<SAMRAI::hier::RefineOperator> const& refineOp,
                         Key const& key,
                         std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
                         std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr);


        /*@brief add a static communication between sources and destinations.
         * This overload takes several sources/destinations/keys and add one refiner for each*/
        template<typename Resources, typename Keys>
        void
        addStaticRefiners(Resources const& destinations, Resources const& sources,
                          std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, Keys const& keys,
                          std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
                          std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr);


        /*@brief convenience overload of the above when source = destination, for VecField*/
        template<typename Srcs, typename Keys>
        void
        addStaticRefiners(Srcs const& src_dest,
                          std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, Keys const& key,
                          std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
                          std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr);


        // this overload takes simple strings.
        void
        addTimeRefiner(std::string const& ghost, std::string const& model,
                       std::string const& oldModel,
                       std::shared_ptr<SAMRAI::hier::RefineOperator> const& refineOp,
                       std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp,
                       std::string const& key,
                       std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
                       std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr);

        /**
         * @brief fill the given pool of refiners with a new refiner per VecField
         * in ghostVecs. Data will be spatially refined using the specified refinement
         * operator, and time interpolated between time n and n+1 of next coarser data,
         * represented by modelVec and oldModelVec.*/
        void
        addTimeRefiners(std::vector<std::string> const& ghostVecs, std::string const& modelVec,
                        std::string const& oldModelVec,
                        std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOp,
                        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>& timeOp,
                        std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr,
                        std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat  = nullptr);



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



    private:
        using Qty = std::string;
        std::map<Qty, Refiner_t> refiners_;
        std::shared_ptr<ResourcesManager> rm_;
    };


} // namespace amr
} // namespace PHARE


namespace PHARE::amr
{


template<typename ResourcesManager, RefinerType Type>
template<typename Resource, typename Key>
void RefinerPool<ResourcesManager, Type>::addStaticRefiner(
    Resource const& dst, Resource const& src,
    std::shared_ptr<SAMRAI::hier::RefineOperator> const& refineOp, Key const& key,
    std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern,
    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat)
{
    auto const [it, success]
        = refiners_.insert({key, Refiner_t(dst, src, rm_, refineOp, fillPattern, patchStrat)});

    if (!success)
        throw std::runtime_error(key + " is already registered");
}


template<typename ResourcesManager, RefinerType Type>
template<typename Resource, typename Key>
void RefinerPool<ResourcesManager, Type>::addStaticRefiner(
    Resource const& src_dst, std::shared_ptr<SAMRAI::hier::RefineOperator> const& refineOp,
    Key const& key, std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern,
    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat)
{
    addStaticRefiner(src_dst, src_dst, refineOp, key, fillPattern, patchStrat);
}


template<typename ResourcesManager, RefinerType Type>
template<typename Resources, typename Keys>
void RefinerPool<ResourcesManager, Type>::addStaticRefiners(
    Resources const& destinations, Resources const& sources,
    std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, Keys const& keys,
    std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern,
    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat)
{
    assert(destinations.size() == sources.size());
    assert(destinations.size() == keys.size());

    for (std::size_t i = 0; i < destinations.size(); ++i)
        addStaticRefiner(destinations[i], sources[i], refineOp, keys[i], fillPattern, patchStrat);
}


template<typename ResourcesManager, RefinerType Type>
template<typename Srcs, typename Keys>
void RefinerPool<ResourcesManager, Type>::addStaticRefiners(
    Srcs const& src_dest, std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp, Keys const& keys,
    std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern,
    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat)
{
    addStaticRefiners(src_dest, src_dest, refineOp, keys, fillPattern, patchStrat);
}



template<typename ResourcesManager, RefinerType Type>
void RefinerPool<ResourcesManager, Type>::addTimeRefiner(
    std::string const& ghost, std::string const& model, std::string const& oldModel,
    std::shared_ptr<SAMRAI::hier::RefineOperator> const& refineOp,
    std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> const& timeOp, std::string const& key,
    std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern,
    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat)
{
    auto const [it, success] = refiners_.insert(
        {key, Refiner_t(ghost, model, oldModel, rm_, refineOp, timeOp, fillPattern, patchStrat)});
    if (!success)
        throw std::runtime_error(key + " is already registered");
}


template<typename ResourcesManager, RefinerType Type>
void RefinerPool<ResourcesManager, Type>::addTimeRefiners(
    std::vector<std::string> const& ghostVecs, std::string const& modelVec,
    std::string const& oldModelVec, std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOp,
    std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>& timeOp,
    std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern,
    std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> patchStrat)
{
    for (auto const& ghostVec : ghostVecs)
        addTimeRefiner(ghostVec, modelVec, oldModelVec, refineOp, timeOp, ghostVec, fillPattern,
                       patchStrat);
}


} // namespace PHARE::amr

#endif
