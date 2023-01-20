#include "patch_level_ghost_fill_pattern.hpp"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/MathUtilities.h"


namespace PHARE::amr
{
/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */

PatchLevelGhostFillPattern::PatchLevelGhostFillPattern()
    : d_max_fill_boxes(0)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

PatchLevelGhostFillPattern::~PatchLevelGhostFillPattern() {}

/*
 *************************************************************************
 *
 * computeFillBoxesAndNeighborhoodSets
 *
 *************************************************************************
 */
void PatchLevelGhostFillPattern::computeFillBoxesAndNeighborhoodSets(
    std::shared_ptr<SAMRAI::hier::BoxLevel>& fill_box_level,
    std::shared_ptr<SAMRAI::hier::Connector>& dst_to_fill,
    const SAMRAI::hier::BoxLevel& dst_box_level, const SAMRAI::hier::IntVector& fill_ghost_width,
    bool data_on_patch_border)
{
    TBOX_ASSERT_OBJDIM_EQUALITY2(dst_box_level, fill_ghost_width);

    fill_box_level.reset(new SAMRAI::hier::BoxLevel(dst_box_level.getRefinementRatio(),
                                                    dst_box_level.getGridGeometry(),
                                                    dst_box_level.getMPI()));

    dst_to_fill.reset(
        new SAMRAI::hier::Connector(dst_box_level, *fill_box_level, fill_ghost_width));

    const SAMRAI::hier::BoxContainer& dst_boxes = dst_box_level.getBoxes();

    const int dst_level_num = dst_box_level.getGridGeometry()->getEquivalentLevelNumber(
        dst_box_level.getRefinementRatio());

    SAMRAI::hier::IntVector dst_to_dst_width(fill_ghost_width);
    if (data_on_patch_border)
    {
        dst_to_dst_width += SAMRAI::hier::IntVector::getOne(fill_ghost_width.getDim());
    }

    const SAMRAI::hier::Connector& dst_to_dst = dst_box_level.findConnector(
        dst_box_level, dst_to_dst_width, SAMRAI::hier::CONNECTOR_IMPLICIT_CREATION_RULE, true);

    /*
     * To get the level border, grow each patch box and remove
     * the level from it.
     */
    SAMRAI::hier::LocalId last_id = dst_box_level.getLastLocalId();
    for (SAMRAI::hier::RealBoxConstIterator ni(dst_boxes.realBegin()); ni != dst_boxes.realEnd();
         ++ni)
    {
        const SAMRAI::hier::Box& dst_box = *ni;
        SAMRAI::hier::BoxContainer fill_boxes(SAMRAI::hier::Box::grow(dst_box, fill_ghost_width));
        SAMRAI::hier::Connector::ConstNeighborhoodIterator nabrs
            = dst_to_dst.find(dst_box.getBoxId());
        for (SAMRAI::hier::Connector::ConstNeighborIterator na = dst_to_dst.begin(nabrs);
             na != dst_to_dst.end(nabrs); ++na)
        {
            // SAMRAI default PatchLevelBorderFillPattern from which this code
            // is based seems to overwrite the border node
            // therefore since we remove 'na' from the fill_boxes,
            // we grow na by 1 in the hope that it
            auto tmp = SAMRAI::hier::Box::grow(
                *na, SAMRAI ::hier::IntVector::getOne(fill_ghost_width.getDim()));
            if (dst_box.getBlockId() == na->getBlockId())
            {
                fill_boxes.removeIntersections(tmp);
            }
            else
            {
                throw std::runtime_error("WE SHOULD NEVER GET HERE");
            }
        }

        if (!fill_boxes.empty())
        {
            d_max_fill_boxes
                = SAMRAI::tbox::MathUtilities<int>::Max(d_max_fill_boxes, fill_boxes.size());
            SAMRAI::hier::Connector::NeighborhoodIterator base_box_itr
                = dst_to_fill->makeEmptyLocalNeighborhood(dst_box.getBoxId());
            for (SAMRAI::hier::BoxContainer::iterator li = fill_boxes.begin();
                 li != fill_boxes.end(); ++li)
            {
                SAMRAI::hier::Box fill_box(*li, ++last_id, dst_box.getOwnerRank());
                TBOX_ASSERT(fill_box.getBlockId() == dst_box.getBlockId());
                fill_box_level->addBoxWithoutUpdate(fill_box);
                dst_to_fill->insertLocalNeighbor(fill_box, base_box_itr);
            }
        }
    }
    fill_box_level->finalize();
}

void PatchLevelGhostFillPattern::computeDestinationFillBoxesOnSourceProc(
    FillSet& dst_fill_boxes_on_src_proc, const SAMRAI::hier::BoxLevel& dst_box_level,
    const SAMRAI::hier::Connector& src_to_dst, const SAMRAI::hier::IntVector& fill_ghost_width)
{
    NULL_USE(dst_box_level);
    NULL_USE(src_to_dst);
    NULL_USE(fill_ghost_width);
    NULL_USE(dst_fill_boxes_on_src_proc);
    if (!needsToCommunicateDestinationFillBoxes())
    {
        TBOX_ERROR("PatchLevelGhostFillPattern cannot compute destination:\n"
                   << "fill boxes on the source processor.\n");
    }
}

bool PatchLevelGhostFillPattern::needsToCommunicateDestinationFillBoxes() const
{
    return true;
}

bool PatchLevelGhostFillPattern::doesSourceLevelCommunicateToDestination() const
{
    return false;
}

bool PatchLevelGhostFillPattern::fillingCoarseFineGhosts() const
{
    return true;
}

bool PatchLevelGhostFillPattern::fillingEnhancedConnectivityOnly() const
{
    return false;
}

int PatchLevelGhostFillPattern::getMaxFillBoxes() const
{
    return d_max_fill_boxes;
}
} // namespace PHARE::amr
