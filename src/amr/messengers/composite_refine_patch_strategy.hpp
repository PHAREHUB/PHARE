#ifndef PHARE_AMR_COMPOSITE_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_COMPOSITE_REFINE_PATCH_STRATEGY_HPP

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <memory>
#include <vector>

namespace PHARE::amr
{

/**
 * @brief A SAMRAI RefinePatchStrategy that forwards every hook to a list of strategies.
 *
 * A RefineSchedule takes a single patch strategy, while each PHARE field patch strategy handles
 * the one patch data id it registered. A schedule filling several items gives each item its
 * strategy through this list. The strategies are called in list order; each must touch only its
 * own patch data, so the order does not change the result.
 *
 * The stencil width is the maximum over the strategies, so each strategy runs with the widest one,
 * e.g. the magnetic strategy (stencil 0) runs with the moments' width 1: SAMRAI may then
 * prepare one more coarse layer for it and apply its boundary condition there, which a
 * schedule of the magnetic field alone would not do.
 */
class CompositeRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy
{
public:
    using strategy_ptr = std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy>;

    explicit CompositeRefinePatchStrategy(std::vector<strategy_ptr>&& strategies)
        : strategies_{std::move(strategies)}
    {
    }

    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch& patch, double const fill_time,
                                       SAMRAI::hier::IntVector const& ghost_width_to_fill) override
    {
        for (auto const& strategy : strategies_)
            strategy->setPhysicalBoundaryConditions(patch, fill_time, ghost_width_to_fill);
    }

    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        auto width = SAMRAI::hier::IntVector::getZero(dim);
        for (auto const& strategy : strategies_)
            width.max(strategy->getRefineOpStencilWidth(dim));
        return width;
    }

    void preprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                          SAMRAI::hier::Box const& fine_box,
                          SAMRAI::hier::IntVector const& ratio) override
    {
        for (auto const& strategy : strategies_)
            strategy->preprocessRefine(fine, coarse, fine_box, ratio);
    }

    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
        for (auto const& strategy : strategies_)
            strategy->postprocessRefine(fine, coarse, fine_box, ratio);
    }

    // RefineSchedule calls only the box versions below, never the single-box hooks above. Each
    // strategy then reaches its own single-box hook through its default box implementation.
    void preprocessRefineBoxes(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                               SAMRAI::hier::BoxContainer const& fine_boxes,
                               SAMRAI::hier::IntVector const& ratio) override
    {
        for (auto const& strategy : strategies_)
            strategy->preprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }

    void postprocessRefineBoxes(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                                SAMRAI::hier::BoxContainer const& fine_boxes,
                                SAMRAI::hier::IntVector const& ratio) override
    {
        for (auto const& strategy : strategies_)
            strategy->postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }

    // a host-device synchronization is needed if any strategy asks for one
    void preprocessRefineLevel(
        SAMRAI::hier::PatchLevel& fine_level, SAMRAI::hier::PatchLevel const& coarse_level,
        SAMRAI::hier::Connector const& coarse_to_fine,
        SAMRAI::hier::Connector const& coarse_to_unfilled,
        std::vector<std::vector<std::shared_ptr<SAMRAI::hier::BoxOverlap>>> const& overlaps,
        SAMRAI::xfer::RefineClasses::Data const** refine_items) override
    {
        bool synchronize = false;
        for (auto const& strategy : strategies_)
        {
            strategy->preprocessRefineLevel(fine_level, coarse_level, coarse_to_fine,
                                            coarse_to_unfilled, overlaps, refine_items);
            synchronize = strategy->needSynchronize() or synchronize;
        }
        setNeedRefineSynchronize(synchronize);
    }

    void postprocessRefineLevel(SAMRAI::hier::PatchLevel& fine_level,
                                SAMRAI::hier::PatchLevel const& coarse_level,
                                SAMRAI::hier::Connector const& coarse_to_fine,
                                SAMRAI::hier::Connector const& coarse_to_unfilled) override
    {
        bool synchronize = false;
        for (auto const& strategy : strategies_)
        {
            strategy->postprocessRefineLevel(fine_level, coarse_level, coarse_to_fine,
                                             coarse_to_unfilled);
            synchronize = strategy->needSynchronize() or synchronize;
        }
        setNeedRefineSynchronize(synchronize);
    }

private:
    std::vector<strategy_ptr> strategies_;
};

} // namespace PHARE::amr

#endif
