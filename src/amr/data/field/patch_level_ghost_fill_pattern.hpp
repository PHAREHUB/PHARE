#ifndef PHARE_PATCH_LEVEL_GHOST_FILL_PATTERN_HPP
#define PHARE_PATCH_LEVEL_GHOST_FILL_PATTERN_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"

namespace PHARE::amr
{
/*!
 * @brief PatchLevelFillPattern implementation for filling at PatchLevel
 * boundaries.
 *
 * For documentation on this interface see @ref PatchLevelFillPattern
 *
 * The fill boxes will consist of the ghost regions lying outside of
 * the level interior--in other words the ghost regions at physical
 * and coarse-fine boundaries.
 */

class PatchLevelGhostFillPattern : public SAMRAI::xfer::PatchLevelFillPattern
{
public:
    /*!
     * @brief Default constructor
     */
    PatchLevelGhostFillPattern();

    /*!
     * @brief Destructor
     */
    virtual ~PatchLevelGhostFillPattern();

    /*!
     * @brief Compute the boxes to be filled and related communication data.
     *
     * The computed fill_box_level will cover the ghost regions around the
     * boxes of dst_box_level at coarse-fine and physical boundaries.
     * The width of those ghost regions will be determined by fill_ghost_width.
     *
     * @param[out] fill_box_level       Output BoxLevel to be filled
     * @param[out] dst_to_fill          Output Connector between
     *                                  dst_box_level and fill_box_level
     * @param[in] dst_box_level         destination level
     * @param[in] fill_ghost_width      Ghost width being filled by refine
     *                                  schedule
     * @param[in] data_on_patch_border  true if there is data living on patch
     *                                  borders
     *
     * @pre dst_box_level.getDim() == fill_ghost_width.getDim()
     */
    void
    computeFillBoxesAndNeighborhoodSets(std::shared_ptr<SAMRAI::hier::BoxLevel>& fill_box_level,
                                        std::shared_ptr<SAMRAI::hier::Connector>& dst_to_fill,
                                        const SAMRAI::hier::BoxLevel& dst_box_level,
                                        const SAMRAI::hier::IntVector& fill_ghost_width,
                                        bool data_on_patch_border);

    /*!
     * @brief Return true to indicate source patch owners cannot compute
     * fill boxes without using communication.
     */
    bool needsToCommunicateDestinationFillBoxes() const;

    /*!
     * @brief Virtual method to compute the destination fill boxes.
     *
     * Since needsToCommunicateDestinationFillBoxes() returns true, this
     * method should never be called.  It is implemented here to satisfy
     * the pure virtual interface from the base class.  An error will result
     * if this is ever called.
     *
     * @pre needsToCommunicateDestinationFillBoxes()
     */
    void computeDestinationFillBoxesOnSourceProc(FillSet& dst_fill_boxes_on_src_proc,
                                                 const SAMRAI::hier::BoxLevel& dst_box_level,
                                                 const SAMRAI::hier::Connector& src_to_dst,
                                                 const SAMRAI::hier::IntVector& fill_ghost_width);

    /*!
     * @brief Tell RefineSchedule not to communicate data directly from source
     * to destination level.
     *
     * RefineSchedule should not attempt to fill as much of the
     * destination level as possible from the source level at the same
     * level of resolution for this fill pattern.
     *
     * @return false
     */
    bool doesSourceLevelCommunicateToDestination() const;

    /*!
     * @brief Return the maximum number of fill boxes.
     *
     * @return maximum number of fill boxes.
     */
    int getMaxFillBoxes() const;

    /*!
     * @brief Returns true because this fill pattern fills coarse-fine ghost
     * data.
     */
    bool fillingCoarseFineGhosts() const;

    /*!
     * @brief Returns false because this fill pattern is not specialized for
     * enhanced connectivity only.
     */
    bool fillingEnhancedConnectivityOnly() const;

private:
    PatchLevelGhostFillPattern(const PatchLevelGhostFillPattern&);            // not implemented
    PatchLevelGhostFillPattern& operator=(const PatchLevelGhostFillPattern&); // not implemented

    /*!
     * @brief Maximum number of fill boxes across all destination patches.
     */
    int d_max_fill_boxes;
};
} // namespace PHARE::amr

#endif
