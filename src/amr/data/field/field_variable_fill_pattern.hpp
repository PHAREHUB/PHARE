#ifndef PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP


#include "core/def/phare_mpi.hpp"
#include <core/hybrid/hybrid_quantities.hpp>

#include <amr/utilities/box/amr_box.hpp>
#include "amr/data/field/field_geometry.hpp"

#include <SAMRAI/pdat/CellOverlap.h>
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <cassert>

namespace PHARE::amr
{
/*
  This class is used from multiple schedules
  To know which schedule we are coming from, we have `std::optional<bool> opt_overwrite_interior_`

  the modes are :

    1. To synchronize primal nodes on shared patch borders
        e.g. hybrid_hybrid_messenger_strategy.hpp
        HybridHybridMessengerStrategy::magneticSharedNodes_

        in this case, the fillPattern is constructed
        with "std::optional<bool> opt_overwrite_interior_ == std::nullopt",
        we set the forwarding flag of "bool overwrite_interior" to true by default
        and it is only set to false for one of the 2 patches involved in the overlap
        so that only one process assigns its value to the shared border node
        We also remove the exclusive interior of the src patch to isolate only shared primal
  nodes.

    2. To synchronize pure ghost values from src domain values
    in that case, the optional is set to "false" and overwrite_interior takes this value
    none of the two patches overwrites the shared border nodes and only pure ghost nodes are
    filled.

  Notes on shared-node overwrite interior: https://github.com/LLNL/SAMRAI/issues/170

*/
// This class is mostly a copy of BoxGeometryVariableFillPattern
template<std::size_t dimension>
class FieldFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    constexpr static std::size_t dim = dimension;

public:
    FieldFillPattern(std::optional<bool> overwrite_interior = false)
        : opt_overwrite_interior_{overwrite_interior}
    {
    }


    virtual ~FieldFillPattern() {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const fn_overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
#ifndef DEBUG_CHECK_DIM_ASSERTIONS
        NULL_USE(dst_patch_box);
#endif
        TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);

        bool overwrite_interior = true; // replace func param
        assert(fn_overwrite_interior == overwrite_interior);

        if (opt_overwrite_interior_) // not node only
        {
            // this sets overwrite_interior to false
            overwrite_interior = *opt_overwrite_interior_;
            assert(overwrite_interior == false); // never true
        }

        // overwrite_interior is always false here
        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior,
                                             transformation);
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    FieldFillPattern(FieldFillPattern const&)            = delete;
    FieldFillPattern& operator=(FieldFillPattern const&) = delete;

    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
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
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& node_fill_boxes,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
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

    std::optional<bool> opt_overwrite_interior_{std::nullopt};
};




template<typename Gridlayout_t> // ASSUMED ALL PRIMAL!
class FieldGhostInterpOverlapFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    std::size_t constexpr static dim = Gridlayout_t::dimension;
    using FieldGeometry_t            = FieldGeometryBase<dim>;

public:
    FieldGhostInterpOverlapFillPattern() {}
    ~FieldGhostInterpOverlapFillPattern() override {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        PHARE_LOG_SCOPE(3, "FieldGhostInterpOverlapFillPattern::calculateOverlap");

        if (phare_box_from<dim>(dst_patch_box) == phare_box_from<dim>(src_mask))
            return std::make_shared<FieldOverlap>(SAMRAI::hier::BoxContainer{}, transformation);

        auto const _primal_ghost_box = [](auto const& box) {
            auto gb = grow(box, Gridlayout_t::nbrGhosts());
            gb.upper += 1;
            return gb;
        };

        auto const src_ghost_box
            = core::shift(_primal_ghost_box(phare_box_from<dim>(
                              dynamic_cast<FieldGeometry_t const&>(src_geometry).patchBox)),
                          as_point<dim>(transformation));

        auto const dst_ghost_box = _primal_ghost_box(
            phare_box_from<dim>(dynamic_cast<FieldGeometry_t const&>(dst_geometry).patchBox));

        SAMRAI::hier::BoxContainer dest;
        if (auto overlap = dst_ghost_box * src_ghost_box)
            dest.push_back(samrai_box_from(*overlap));

        return std::make_shared<FieldOverlap>(dest, transformation);
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        throw std::runtime_error("never called");
    }


    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& node_fill_boxes,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    {
        throw std::runtime_error("no refinement supported or expected");
    }
};



} // namespace PHARE::amr

#endif /* PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_H */
