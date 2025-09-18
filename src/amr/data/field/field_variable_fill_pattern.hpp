#ifndef PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP

#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp"
#include "core/utilities/types.hpp"
#include <core/hybrid/hybrid_quantities.hpp>
#include "core/data/tensorfield/tensorfield.hpp"

#include <amr/utilities/box/amr_box.hpp>
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"
#include "amr/data/tensorfield/tensor_field_geometry.hpp"

#include <SAMRAI/pdat/CellOverlap.h>
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <cassert>
#include <cstddef>

namespace PHARE::amr
{
/*
  This class is used to synchronize pure ghost values from src domain values
    in that case, we default overwrite_interior to "false" so none of the two patches overwrites
    the shared border nodes and only pure ghost nodes are filled.

  Notes on shared-node overwrite interior: https://github.com/LLNL/SAMRAI/issues/170
*/
// This class is mostly a copy of BoxGeometryVariableFillPattern
template<std::size_t dimension>
class FieldFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    constexpr static std::size_t dim = dimension;

public:
    // defaulted param makes this the default constructor also
    FieldFillPattern(bool overwrite_interior = false)
        : overwrite_interior_{overwrite_interior}

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

        assert(fn_overwrite_interior == true); // expect default as true

        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior_,

                                             transformation);
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


        SAMRAI::hier::Transformation transformation(
            SAMRAI::hier::IntVector::getZero(patch_box.getDim()));

        SAMRAI::hier::BoxContainer overlap_boxes(fill_boxes);
        overlap_boxes.intersectBoxes(data_box);

        auto geom = pdf.getBoxGeometry(patch_box);
        auto basic_overlap
            = pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes, transformation);

        if (overwrite_interior_)
            return basic_overlap;

        // from here we do not overwrite interior values
        // so we need to remove from the overlap boxes their intersection
        // with the interior box.
        // Note this only removes cells that belong to A PATCH interior
        // i.e. some cells in this overlap may still be level interior cell but
        // belong to another patch interior.
        auto& overlap         = dynamic_cast<FieldOverlap const&>(*basic_overlap);
        auto destinationBoxes = overlap.getDestinationBoxContainer();
        auto& casted          = dynamic_cast<FieldGeometryBase<dimension> const&>(*geom);
        destinationBoxes.removeIntersections(casted.interiorFieldBox());

        return std::make_shared<FieldOverlap>(destinationBoxes, overlap.getTransformation());
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

    bool const overwrite_interior_;
};


template<std::size_t dimension, std::size_t rank_ = 1>
class TensorFieldFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    static constexpr std::size_t N = core::detail::tensor_field_dim_from_rank<rank_>();

public:
    TensorFieldFillPattern(bool overwrite_interior = false)
        : overwrite_interior_{overwrite_interior}
    {
    }

    ~TensorFieldFillPattern() override = default;

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(const SAMRAI::hier::BoxGeometry& dst_geometry,
                     const SAMRAI::hier::BoxGeometry& src_geometry,
                     const SAMRAI::hier::Box& dst_patch_box, const SAMRAI::hier::Box& src_mask,
                     const SAMRAI::hier::Box& fill_box, bool const fn_overwrite_interior,
                     const SAMRAI::hier::Transformation& transformation) const override
    {
        // Note fn_overwrite_interior is the boolean passed by SAMRAI and is always true
        // this `VariableFillPattern` overrides this behavior using its own `overwrite_interior_`
        // set on construction and depending on the use case.
        return dst_geometry.calculateOverlap(src_geometry, src_mask, fill_box, overwrite_interior_,
                                             transformation);
    }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& fill_boxes,
                            SAMRAI::hier::BoxContainer const& node_fill_boxes,
                            SAMRAI::hier::Box const& patch_box, SAMRAI::hier::Box const& data_box,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    {
        SAMRAI::hier::Transformation transformation(
            SAMRAI::hier::IntVector::getZero(patch_box.getDim()));

        SAMRAI::hier::BoxContainer overlap_boxes(fill_boxes);
        overlap_boxes.intersectBoxes(data_box);

        auto basic_overlap
            = pdf.getBoxGeometry(patch_box)->setUpOverlap(overlap_boxes, transformation);

        if (overwrite_interior_)
            return basic_overlap;

        // from here we do not overwrite interior values
        // so we need to remove from the overlap boxes their intersection
        // with the interior box.
        // Note this only removes cells that belong to A PATCH interior
        // i.e. some cells in this overlap may still be level interior cell but
        // belong to another patch interior.
        auto geom      = pdf.getBoxGeometry(patch_box);
        auto& casted   = dynamic_cast<TensorFieldGeometryBase<dimension> const&>(*geom);
        auto& toverlap = dynamic_cast<TensorFieldOverlap<rank_> const&>(*basic_overlap);
        auto&& interiorTensorFieldBox = casted.interiorTensorFieldBox();

        auto overlaps = core::for_N<N, core::for_N_R_mode::make_array>([&](auto i) {
            auto& overlap          = toverlap[i];
            auto& interiorFieldBox = interiorTensorFieldBox[i];
            auto destinationBoxes  = overlap->getDestinationBoxContainer();
            destinationBoxes.removeIntersections(interiorFieldBox);

            return std::make_shared<FieldOverlap>(destinationBoxes, overlap->getTransformation());
        });

        return std::make_shared<TensorFieldOverlap<rank_>>(std::move(overlaps));
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    TensorFieldFillPattern(TensorFieldFillPattern const&)            = delete;
    TensorFieldFillPattern& operator=(TensorFieldFillPattern const&) = delete;

    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        TBOX_ERROR("getStencilWidth() should not be called for TensorFieldFillPattern.");
        return SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension(1));
    }

    bool overwrite_interior_;
};


// We use this fill pattern to sum the contributions of border fields like rho and flux
/** \brief VariableFillPattern that is used to fill incomplete ghost domain moment nodes
 *
 * After deposition of domain particles, some domain and ghost nodes lack contributions
 * from particle that exist on a neighboring patch.
 * The extent of incomplete nodes in the ghost layer and in domain depends on interp order.
 *
 * Example, at interpolation order 1, only the border node will be incomplete after
 * depositing domain particles since these hit only the two primal nodes surronding its position.
 * However, we deposit also leaving domain particles before they are sent to patchGhost particles
 * and shipped to neighrboring patches.
 * Leaving particles can be found in the first ghost cell from domain, so the first primal
 * ghost node from domain will also have some incomplete contribution.
 *
 * At order 1, thus, there is an overlap of 3 primal nodes that are incomplete:
 * the shared border node and the first domain and first ghost.
 *
 *                        ghost cells <-|-> patch
 *                          +           +           +
 *                          |  leaving  | domain particles
 *                          | particles |
 *
 *
 * As a first try and to keep it simple, this fill pattern simply creates the overlap
 * that is the intersection of the field ghost boxes of the source and destination patch datas.
 * That is, at interpolation 1 we have 2 ghost cells  thus it is 5 nodes that overlap
 * even though the outermost ghost should have 0 contribution.
 *
 *                        ghost cells <-|-> patch
 *               +          +           +           +          +
 *               ^          |  leaving  | domain particles
 *               |          | particles |
 *               0
 * */
template<typename Gridlayout_t> // ASSUMED ALL PRIMAL!
class FieldGhostInterpOverlapFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    std::size_t constexpr static dim = Gridlayout_t::dimension;
    using FieldGeometry_t            = FieldGeometryBase<dim>;

public:
    FieldGhostInterpOverlapFillPattern() {}
    ~FieldGhostInterpOverlapFillPattern() override {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& _dst_geometry,
                     SAMRAI::hier::BoxGeometry const& _src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        PHARE_LOG_SCOPE(3, "FieldGhostInterpOverlapFillPattern::calculateOverlap");

        // Skip if src and dst are the same
        if (phare_box_from<dim>(dst_patch_box) == phare_box_from<dim>(src_mask))
            return std::make_shared<FieldOverlap>(SAMRAI::hier::BoxContainer{}, transformation);

        if (dynamic_cast<FieldGeometry_t const*>(&_dst_geometry))
            return calculateOverlap(dynamic_cast<FieldGeometry_t const&>(_dst_geometry),
                                    dynamic_cast<FieldGeometry_t const&>(_src_geometry),
                                    dst_patch_box, src_mask, fill_box, overwrite_interior,
                                    transformation);
        else
            throw std::runtime_error("bad cast");
    }


    std::shared_ptr<SAMRAI::hier::BoxOverlap> static calculateOverlap(
        auto const& dst_geometry, auto const& src_geometry, SAMRAI::hier::Box const& dst_patch_box,
        SAMRAI::hier::Box const& src_mask, SAMRAI::hier::Box const& fill_box,
        bool const overwrite_interior, SAMRAI::hier::Transformation const& transformation)
    {
        auto const _primal_ghost_box = [](auto const& box) {
            auto gb = grow(box, Gridlayout_t::nbrGhosts());
            gb.upper += 1;
            return gb;
        };

        auto const src_ghost_box = [&]() {
            auto const box              = phare_box_from<dim>(src_geometry.patchBox);
            auto const primal_ghost_box = _primal_ghost_box(box);
            return amr::shift(primal_ghost_box, transformation);
        }();

        auto const dst_ghost_box = [&]() {
            auto const box = phare_box_from<dim>(dst_geometry.patchBox);
            return _primal_ghost_box(box);
        }();


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

template<typename Gridlayout_t, std::size_t rank_ = 1> // ASSUMED ALL PRIMAL!
class TensorFieldGhostInterpOverlapFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    std::size_t constexpr static dim = Gridlayout_t::dimension;
    static constexpr auto N          = core::detail::tensor_field_dim_from_rank<rank_>();

    using TensorFieldGeometry_t = TensorFieldGeometryBase<dim, rank_>;

public:
    TensorFieldGhostInterpOverlapFillPattern() {}
    ~TensorFieldGhostInterpOverlapFillPattern() override {}

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& _dst_geometry,
                     SAMRAI::hier::BoxGeometry const& _src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        PHARE_LOG_SCOPE(3, "TensorFieldGhostInterpOverlapFillPattern::calculateOverlap");

        // Skip if src and dst are the same
        if (phare_box_from<dim>(dst_patch_box) == phare_box_from<dim>(src_mask))
        {
            auto overlaps = core::for_N<N, core::for_N_R_mode::make_array>([&](auto /*i*/) {
                return std::make_shared<FieldOverlap>(SAMRAI::hier::BoxContainer{}, transformation);
            });
            return std::make_shared<TensorFieldOverlap<rank_>>(std::move(overlaps));
        }

        if (dynamic_cast<TensorFieldGeometry_t const*>(&_dst_geometry))
        {
            auto overlaps = core::for_N<N, core::for_N_R_mode::make_array>([&](auto /*i*/) {
                auto overlap = FieldGhostInterpOverlapFillPattern<Gridlayout_t>::calculateOverlap(
                    dynamic_cast<TensorFieldGeometry_t const&>(_dst_geometry),
                    dynamic_cast<TensorFieldGeometry_t const&>(_src_geometry), dst_patch_box,
                    src_mask, fill_box, overwrite_interior, transformation);

                return std::dynamic_pointer_cast<FieldOverlap>(overlap);
            });
            return std::make_shared<TensorFieldOverlap<rank_>>(std::move(overlaps));
        }

        else
            throw std::runtime_error("bad cast");
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
