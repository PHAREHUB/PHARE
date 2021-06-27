
#ifndef GhostOnlyVariablyFillPattern_
#define GhostOnlyVariablyFillPattern_

#include <cassert>
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "amr/data/field/field_geometry.h"

namespace PHARE::amr
{
class GhostOnlyVariablyFillPattern : public SAMRAI::xfer::VariableFillPattern
{
public:
    GhostOnlyVariablyFillPattern() {}
    virtual ~GhostOnlyVariablyFillPattern() {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(const SAMRAI::hier::BoxGeometry& dst_geometry,
                     const SAMRAI::hier::BoxGeometry& src_geometry,
                     const SAMRAI::hier::Box& dst_patch_box, const SAMRAI::hier::Box& src_mask,
                     const SAMRAI::hier::Box& fill_box, const bool overwrite_interior_,
                     const SAMRAI::hier::Transformation& transformation) const
    {
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
        NULL_USE(dst_patch_box);
#endif
        TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);

        bool overwrite_interior = true; // replace func param
        assert(overwrite_interior_ == overwrite_interior);

        auto& dst_cast = dynamic_cast<AFieldGeometry const&>(dst_geometry);
        auto& src_cast = dynamic_cast<AFieldGeometry const&>(src_geometry);

        // for shared border node value sync
        if (src_cast.patchBox.getGlobalId().getOwnerRank()
            != dst_cast.patchBox.getGlobalId().getOwnerRank())
            overwrite_interior = src_cast.patchBox.getLocalId() > dst_cast.patchBox.getLocalId();

        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior,
                                             transformation);
    }

    const std::string& getPatternName() const { return s_name_id; }

private:
    GhostOnlyVariablyFillPattern(const GhostOnlyVariablyFillPattern&);            // not implemented
    GhostOnlyVariablyFillPattern& operator=(const GhostOnlyVariablyFillPattern&); // not implemented

    static const std::string s_name_id; // = "GHOST_ONLY_FILL_PATTERN";

    const SAMRAI::hier::IntVector& getStencilWidth()
    {
        TBOX_ERROR("getStencilWidth() should not be\n"
                   << "called.  This pattern creates overlaps based on\n"
                   << "the BoxGeometry objects and is not restricted to a\n"
                   << "specific stencil.\n");

        /*
         * Dummy return value that will never get reached.
         */
        return SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(1));
    }

    /*
     *************************************************************************
     *
     * Compute BoxOverlap that specifies data to be filled by refinement
     * operator.
     *
     *************************************************************************
     */
    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(const SAMRAI::hier::BoxContainer& fill_boxes,
                            const SAMRAI::hier::BoxContainer& node_fill_boxes,
                            const SAMRAI::hier::Box& patch_box, const SAMRAI::hier::Box& data_box,
                            const SAMRAI::hier::PatchDataFactory& pdf) const
    {
        NULL_USE(node_fill_boxes);

        /*
         * For this (default) case, the overlap is simply the intersection of
         * fill_boxes and data_box.
         */
        SAMRAI::hier::Transformation transformation(
            SAMRAI::hier::IntVector::getZero(patch_box.getDim()));

        SAMRAI::hier::BoxContainer overlap_boxes(fill_boxes);
        overlap_boxes.intersectBoxes(data_box);

        return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes, transformation);
    }
};

} // namespace PHARE::amr

#endif /*GhostOnlyVariablyFillPattern*/
