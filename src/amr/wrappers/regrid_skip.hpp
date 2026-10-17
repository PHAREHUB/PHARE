#ifndef PHARE_AMR_WRAPPERS_REGRID_SKIP_HPP
#define PHARE_AMR_WRAPPERS_REGRID_SKIP_HPP

#include "core/def.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep

#include <SAMRAI/hier/BoxLevel.h>
#include <SAMRAI/hier/Connector.h>
#include <SAMRAI/hier/OverlapConnectorAlgorithm.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/mesh/BoxGeneratorStrategy.h>
#include <SAMRAI/mesh/LoadBalanceStrategy.h>

#include <memory>
#include <vector>


namespace PHARE::amr
{
/*
 * SAMRAI recreates the finest level at every regrid point, which with refined timestepping
 * happens every 2 steps of the next coarser level, and is not configurable.
 * Every other regrid of the finest level, we give SAMRAI back the boxes it already has
 * and prevent the load balancer from changing them, such that the new level is identical
 * to the old one and its patch datas can be moved over, see
 * MultiPhysicsIntegrator::initializeLevelData
 */
struct RegridSkipState
{
    int skip_ln = -1;                   // level being recreated with unchanged boxes
    std::vector<std::size_t> regrids{}; // per level count of finest level regrids
};


class RegridSkipBoxGenerator : public SAMRAI::mesh::BoxGeneratorStrategy
{
public:
    RegridSkipBoxGenerator(std::shared_ptr<SAMRAI::mesh::BoxGeneratorStrategy> generator,
                           std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
                           std::shared_ptr<RegridSkipState> state)
        : generator_{std::move(generator)}
        , hierarchy_{std::move(hierarchy)}
        , state_{std::move(state)}
    {
        state_->regrids.resize(hierarchy_->getMaxNumberOfLevels());
    }

    void findBoxesContainingTags(std::shared_ptr<SAMRAI::hier::BoxLevel>& new_box_level,
                                 std::shared_ptr<SAMRAI::hier::Connector>& tag_to_new,
                                 std::shared_ptr<SAMRAI::hier::PatchLevel> const& tag_level,
                                 int const tag_data_index, int const tag_val,
                                 SAMRAI::hier::BoxContainer const& bound_boxes,
                                 SAMRAI::hier::IntVector const& min_box,
                                 SAMRAI::hier::IntVector const& tag_to_new_width) override
    {
        auto const new_ln = tag_level->getLevelNumber() + 1;
        state_->skip_ln   = should_skip_(new_ln) ? new_ln : -1;

        if (state_->skip_ln == new_ln)
            reuse_boxes_(new_box_level, tag_to_new, *tag_level, tag_to_new_width);
        else
            generator_->findBoxesContainingTags(new_box_level, tag_to_new, tag_level,
                                                tag_data_index, tag_val, bound_boxes, min_box,
                                                tag_to_new_width);
    }

    void setMinimumCellRequest(std::size_t const num_cells) override
    {
        generator_->setMinimumCellRequest(num_cells);
    }

    void setRatioToNewLevel(SAMRAI::hier::IntVector const& ratio) override
    {
        generator_->setRatioToNewLevel(ratio);
    }

private:
    // only the finest level, the coarser ones may need to grow to nest finer levels
    bool should_skip_(int const new_ln)
    {
        bool const is_finest_regrid
            = hierarchy_->levelExists(new_ln) and !hierarchy_->levelExists(new_ln + 1);
        return is_finest_regrid and state_->regrids[new_ln]++ % 2 == 1;
    }

    // new level boxes are the current ones, in tag level index space, with the same owners
    void reuse_boxes_(std::shared_ptr<SAMRAI::hier::BoxLevel>& new_box_level,
                      std::shared_ptr<SAMRAI::hier::Connector>& tag_to_new,
                      SAMRAI::hier::PatchLevel const& tag_level,
                      SAMRAI::hier::IntVector const& tag_to_new_width) const
    {
        auto const new_ln         = tag_level.getLevelNumber() + 1;
        auto const& tag_box_level = *tag_level.getBoxLevel();
        auto const& ratio         = hierarchy_->getRatioToCoarserLevel(new_ln);

        new_box_level = std::make_shared<SAMRAI::hier::BoxLevel>(tag_box_level.getRefinementRatio(),
                                                                 tag_box_level.getGridGeometry(),
                                                                 tag_box_level.getMPI());

        for (auto box : hierarchy_->getBoxLevel(new_ln)->getBoxes())
            if (!box.isPeriodicImage())
            {
                box.coarsen(ratio);
                new_box_level->addBoxWithoutUpdate(box);
            }
        new_box_level->finalize();

        // periodic relationships are added by the GriddingAlgorithm later
        SAMRAI::hier::OverlapConnectorAlgorithm{}.findOverlapsWithTranspose(
            tag_to_new, tag_box_level, *new_box_level, tag_to_new_width, tag_to_new_width);
        tag_to_new->removePeriodicRelationships();
        tag_to_new->getTranspose().removePeriodicRelationships();
    }

    std::shared_ptr<SAMRAI::mesh::BoxGeneratorStrategy> generator_;
    std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy_;
    std::shared_ptr<RegridSkipState> state_;
};


class RegridSkipLoadBalancer : public SAMRAI::mesh::LoadBalanceStrategy
{
public:
    RegridSkipLoadBalancer(std::shared_ptr<SAMRAI::mesh::LoadBalanceStrategy> balancer,
                           std::shared_ptr<RegridSkipState> state)
        : balancer_{std::move(balancer)}
        , state_{std::move(state)}
    {
    }

    NO_DISCARD bool getLoadBalanceDependsOnPatchData(int const level_number) const override
    {
        return balancer_->getLoadBalanceDependsOnPatchData(level_number);
    }

    void loadBalanceBoxLevel(SAMRAI::hier::BoxLevel& balance_box_level,
                             SAMRAI::hier::Connector* balance_to_anchor,
                             std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                             int const level_number, SAMRAI::hier::IntVector const& min_size,
                             SAMRAI::hier::IntVector const& max_size,
                             SAMRAI::hier::BoxLevel const& domain_box_level,
                             SAMRAI::hier::IntVector const& bad_interval,
                             SAMRAI::hier::IntVector const& cut_factor,
                             SAMRAI::tbox::RankGroup const& rank_group) const override
    {
        if (level_number == state_->skip_ln)
            return; // keep boxes and owners

        balancer_->loadBalanceBoxLevel(balance_box_level, balance_to_anchor, hierarchy,
                                       level_number, min_size, max_size, domain_box_level,
                                       bad_interval, cut_factor, rank_group);
    }

    void setWorkloadPatchDataIndex(int const data_id, int const level_number) override
    {
        balancer_->setWorkloadPatchDataIndex(data_id, level_number);
    }

private:
    std::shared_ptr<SAMRAI::mesh::LoadBalanceStrategy> balancer_;
    std::shared_ptr<RegridSkipState> state_;
};

} // namespace PHARE::amr

#endif /* PHARE_AMR_WRAPPERS_REGRID_SKIP_HPP */
