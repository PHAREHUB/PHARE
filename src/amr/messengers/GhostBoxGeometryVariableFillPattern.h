
#ifndef GhostBoxGeometryVariableFillPattern_
#define GhostBoxGeometryVariableFillPattern_

#include "SAMRAI/xfer/VariableFillPattern.h"

namespace SAMRAI::xfer
{
class GhostBoxGeometryVariableFillPattern : public VariableFillPattern
{
public:
    /*!
     * @brief Default constructor
     */
    GhostBoxGeometryVariableFillPattern() {}

    /*!
     * @brief Destructor
     */
    virtual ~GhostBoxGeometryVariableFillPattern() {}

    /*!
     * @brief Calculate overlap between the destination and source geometries
     * using the geometries' own overlap calculation methods.
     *
     * The intersection between the given dst_geometry and src_geometry
     * will be calculated according to the properties of those geometries.
     *
     * @param[in] dst_geometry    geometry object for destination box
     * @param[in] src_geometry    geometry object for source box
     * @param[in] dst_patch_box   box for the destination patch
     * @param[in] src_mask        the source mask, the box resulting from
     *                            transforming the source box
     * @param[in] fill_box        the box to be filled
     * @param[in] overwrite_interior  controls whether or not to include the
     *                                destination box interior in the overlap.
     * @param[in] transformation  the transformation from source to
     *                            destination index space.
     *
     * @return                    std::shared_ptr to the calculated overlap
     *                            object
     *
     * @pre dst_patch_box.getDim() == src_mask.getDim()
     */
    std::shared_ptr<hier::BoxOverlap>
    calculateOverlap(const hier::BoxGeometry& dst_geometry, const hier::BoxGeometry& src_geometry,
                     const hier::Box& dst_patch_box, const hier::Box& src_mask,
                     const hier::Box& fill_box, const bool overwrite_interior,
                     const hier::Transformation& transformation) const
    {
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
        NULL_USE(dst_patch_box);
#endif
        TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);
        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box,
                                             /*overwrite_interior*/ false, transformation);
    }

    /*!
     * Computes a BoxOverlap object which defines the space to be filled by
     * a refinement operation.  For this implementation, that space is the
     * intersection between fill_boxes (computed by the RefineSchedule) and
     * data_box, which specifies the extent of the destination data.  The
     * patch data factory is used to compute the overlap with the appropriate
     * data centering, consistent with the centering of the data to be filled.
     *
     * @param[in] fill_boxes  list representing the all of the space on a patch
     *                        or its ghost region that may be filled by a
     *                        refine operator (cell-centered representation)
     * @param[in] node_fill_boxes node-centered representation of fill_boxes
     * @param[in] patch_box   box representing the patch where a refine operator
     *                        will fill data.  (cell-centered representation)
     * @param[in] data_box    box representing the full extent of the region
     *                        covered by a patch data object, including all
     *                        ghosts (cell-centered representation)
     * @param[in] pdf         patch data factory for the data that is to be
     *                        filled
     */
    // std::shared_ptr<hier::BoxOverlap>
    // computeFillBoxesOverlap(const hier::BoxContainer& fill_boxes,
    //                         const hier::BoxContainer& node_fill_boxes, const hier::Box&
    //                         patch_box, const hier::Box& data_box, const hier::PatchDataFactory&
    //                         pdf) const;

    /*!
     * @brief Implementation of interface to get stencil width of a
     * VariableFillPattern.
     *
     * For this class GhostBoxGeometryVariableFillPattern, this method should
     * never be called, since overlaps are computed based on BoxGeometry
     * objects and not on any stencil.  An error will result if this method
     * is invoked.
     */
    // const hier::IntVector& getStencilWidth();

    /*!
     * @brief Returns a string name identifier "BOX_GEOMETRY_FILL_PATTERN".
     */
    const std::string& getPatternName() const { return s_name_id; }

private:
    GhostBoxGeometryVariableFillPattern(
        const GhostBoxGeometryVariableFillPattern&); // not implemented
    GhostBoxGeometryVariableFillPattern&
    operator=(const GhostBoxGeometryVariableFillPattern&); // not implemented

    /*!
     * Static string containing string name identifier for this class
     */
    static const std::string s_name_id; // = "BOX_GEOMETRY_FILL_PATTERN";




    // GhostBoxGeometryVariableFillPattern() {}

    /*
     *************************************************************************
     *
     * Destructor
     *
     *************************************************************************
     */

    // ~GhostBoxGeometryVariableFillPattern() {}

    /*
     *************************************************************************
     *
     * getStencilWidth() throws an error if called.  Only overridding
     * versions of this method in concrete subclasses should be called.
     *
     *************************************************************************
     */
    const hier::IntVector& getStencilWidth()
    {
        TBOX_ERROR("getStencilWidth() should not be\n"
                   << "called.  This pattern creates overlaps based on\n"
                   << "the BoxGeometry objects and is not restricted to a\n"
                   << "specific stencil.\n");

        /*
         * Dummy return value that will never get reached.
         */
        return hier::IntVector::getZero(tbox::Dimension(1));
    }

    /*
     *************************************************************************
     *
     * Compute BoxOverlap that specifies data to be filled by refinement
     * operator.
     *
     *************************************************************************
     */
    std::shared_ptr<hier::BoxOverlap>
    computeFillBoxesOverlap(const hier::BoxContainer& fill_boxes,
                            const hier::BoxContainer& node_fill_boxes, const hier::Box& patch_box,
                            const hier::Box& data_box, const hier::PatchDataFactory& pdf) const
    {
        NULL_USE(node_fill_boxes);

        /*
         * For this (default) case, the overlap is simply the intersection of
         * fill_boxes and data_box.
         */
        hier::Transformation transformation(hier::IntVector::getZero(patch_box.getDim()));

        hier::BoxContainer overlap_boxes(fill_boxes);
        overlap_boxes.intersectBoxes(data_box);

        return pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes, transformation);
    }
};

} // namespace SAMRAI::xfer

#endif /*GhostBoxGeometryVariableFillPattern*/
