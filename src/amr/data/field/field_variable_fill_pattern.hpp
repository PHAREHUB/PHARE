#ifndef PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP
#define PHARE_SRC_AMR_FIELD_FIELD_VARIABLE_FILL_PATTERN_HPP

#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "core/utilities/types.hpp"

#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"
#include "amr/data/tensorfield/tensor_field_geometry.hpp"

#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/hier/BaseGridGeometry.h>
#include <SAMRAI/hier/BoxLevel.h>
#include <SAMRAI/hier/PatchLevel.h>
#include "SAMRAI/xfer/VariableFillPattern.h"

#include <cassert>
#include <cstddef>

namespace PHARE::amr
{
/*
  Mostly a copy of BoxGeometryVariableFillPattern. With overwrite_interior false (the default),
  it fills pure ghost values from src domain values: neither patch overwrites the shared border
  nodes. Notes on shared-node overwrite interior: https://github.com/LLNL/SAMRAI/issues/170

  overwrite_interior is a template parameter so that the two settings are two types, which SAMRAI
  never puts in one equivalence class (it computes the overlaps of a class from its first item).
  Splitting by fill pattern type is fine for this flag only: on the coarse interpolation level,
  where SAMRAI replaces the fill patterns (centered_factory.hpp), the whole patch data is filled
  and the flag does not matter.
*/
template<std::size_t dimension, bool overwrite_interior = false>
class FieldFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    constexpr static std::size_t dim = dimension;

public:
    FieldFillPattern() = default;

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

        return dynamic_cast<FillPatternGeometry const&>(dst_geometry)
            .calculateUnroundedOverlap(src_geometry, src_mask, fill_box, overwrite_interior,
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

        if constexpr (overwrite_interior)
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
};


// overwrite_interior as in FieldFillPattern
template<std::size_t dimension, std::size_t rank_ = 1, bool overwrite_interior = false>
class TensorFieldFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    static constexpr std::size_t N = core::detail::tensor_field_dim_from_rank<rank_>();

public:
    TensorFieldFillPattern() = default;

    ~TensorFieldFillPattern() override = default;

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& dst_patch_box, SAMRAI::hier::Box const& src_mask,
                     SAMRAI::hier::Box const& fill_box, bool const fn_overwrite_interior,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        // Note fn_overwrite_interior is the boolean passed by SAMRAI and is always true
        // this `VariableFillPattern` overrides this behavior using its own `overwrite_interior`
        // set by its type and depending on the use case.
        return dynamic_cast<FillPatternGeometry const&>(dst_geometry)
            .calculateUnroundedOverlap(src_geometry, src_mask, fill_box, overwrite_interior,
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

        if constexpr (overwrite_interior)
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

        auto overlaps = core::for_N_make_array<N>([&](auto i) {
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
};

template<std::size_t dimension>
class NoOverlapFillPattern : public SAMRAI::xfer::VariableFillPattern
{
public:
    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& /*src_geometry*/,
                     SAMRAI::hier::Box const& /*dst_patch_box*/,
                     SAMRAI::hier::Box const& /*src_mask*/, SAMRAI::hier::Box const& /*fill_box*/,
                     bool const /*overwrite_interior*/,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        return dst_geometry.setUpOverlap(SAMRAI::hier::BoxContainer{}, transformation);
    }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& /*fill_boxes*/,
                            SAMRAI::hier::BoxContainer const& /*node_fill_boxes*/,
                            SAMRAI::hier::Box const& patch_box,
                            SAMRAI::hier::Box const& /*data_box*/,
                            SAMRAI::hier::PatchDataFactory const& pdf) const override
    {
        return pdf.getBoxGeometry(patch_box)->setUpOverlap(
            SAMRAI::hier::BoxContainer{},
            SAMRAI::hier::Transformation{SAMRAI::hier::IntVector::getZero(patch_box.getDim())});
    }

    std::string const& getPatternName() const override { return s_name_id; }

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        return SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension{dimension});
    }

private:
    static inline std::string const s_name_id = "NO_OVERLAP_FILL_PATTERN";
};

/** @brief True if both geometries are built on the same patch, not on a periodic image of it.
 * Distinct overlapping patches may have identical boxes, so box ids are compared, not extents.
 */
template<typename Geometry_t>
NO_DISCARD bool isSamePatch(Geometry_t const& dst_geometry, Geometry_t const& src_geometry)
{
    auto const& dst_id = dst_geometry.patchBox.getBoxId();
    auto const& src_id = src_geometry.patchBox.getBoxId();
    if (!dst_id.getLocalId().isValid() or !src_id.getLocalId().isValid())
        throw std::runtime_error("isSamePatch: geometry box without id");
    return dst_id == src_id;
}


/** @brief The boxes of a level, with their periodic images, for OwnerFillPattern.
 *
 * A node held by several patches is owned by the patch with the smallest GlobalId, the rule of
 * makeNonOwnedBoxesFor (amr_utils.hpp). A patch whose box spans a periodic direction holds its
 * node twice, at index i and i + period: the copy with the smallest index owns it. Built on all
 * ranks at once (the boxes are globalized), so all ranks agree without communication.
 */
template<std::size_t dimension>
class SameLevelOwners
{
    using dir_t = SAMRAI::tbox::Dimension::dir_t;

public:
    struct Candidate
    {
        SAMRAI::hier::Box box; // cells, shifted to the periodic image
        SAMRAI::hier::GlobalId id;
        SAMRAI::hier::IntVector shift;
    };

    SameLevelOwners(SAMRAI::hier::PatchLevel const& level,
                    SAMRAI::hier::BaseGridGeometry const& geometry)
    {
        auto const period = geometry.getPeriodicShift(level.getRatioToLevelZero());
        SAMRAI::tbox::Dimension const dim{dimension};

        std::vector<SAMRAI::hier::IntVector> shifts{SAMRAI::hier::IntVector::getZero(dim)};
        for (dir_t d = 0; d < dim.getValue(); ++d)
            if (period(d) != 0)
            {
                auto const unshifted = shifts;
                for (int const sign : {-1, 1})
                    for (auto shift : unshifted)
                    {
                        shift(d) = sign * period(d);
                        shifts.push_back(shift);
                    }
            }

        for (auto const& box : level.getBoxLevel()->getGlobalizedVersion().getGlobalBoxes())
            if (!box.isPeriodicImage())
                for (auto const& shift : shifts)
                {
                    SAMRAI::hier::Box image{box};
                    image.shift(shift);
                    candidates.push_back({image, box.getGlobalId(), shift});
                }
    }

    //! true if `candidate` owns a node before the copy of patch `id` seen through `shift`
    static bool before(Candidate const& candidate, SAMRAI::hier::GlobalId const& id,
                       SAMRAI::hier::IntVector const& shift)
    {
        if (candidate.id != id)
            return candidate.id < id;
        for (dir_t d = 0; d < candidate.shift.getDim().getValue(); ++d)
            if (candidate.shift(d) != shift(d))
                return candidate.shift(d) > shift(d); // its node index is smaller
        return false;
    }

    std::vector<Candidate> candidates;
};


/** @brief Same-level fill where every copy of a node takes the value of its owner's copy.
 *
 * With the default FieldFillPattern, a node shared by patches keeps a copy per patch, and a ghost
 * node of A on the border of B and C is written by both B and C: the last transaction wins. Two
 * copies that differ by round-off then feed asymmetric updates, and the difference grows.
 * Here a node takes the value of its owner (SameLevelOwners) on every patch holding it, the
 * owner's own border nodes and the receiver's included, so no copy may differ after the fill.
 * A transaction only reads owner copies and only writes the other copies: the result does not
 * depend on the order of the transactions.
 *
 * band = false: the nodes of the patches of the level (their border nodes included).
 * band = true: the level ghost band, the ghost nodes outside every patch of the level. A patch
 * holds there the deposit of its level ghost particles, complete one cell outside its box only
 * (hybrid_options.hpp): the owner is the first patch with the node in its box grown by one cell.
 *
 * Same-level schedules only: set the boxes of the level (setOwners) before building one.
 */
template<std::size_t dimension, std::size_t rank, bool band = false>
class OwnerFillPattern : public SAMRAI::xfer::VariableFillPattern
{
    using Owners_t = SameLevelOwners<dimension>;

public:
    void setOwners(std::shared_ptr<Owners_t const> owners) { owners_ = std::move(owners); }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    calculateOverlap(SAMRAI::hier::BoxGeometry const& dst_geometry,
                     SAMRAI::hier::BoxGeometry const& src_geometry,
                     SAMRAI::hier::Box const& /*dst_patch_box*/,
                     SAMRAI::hier::Box const& /*src_mask*/, SAMRAI::hier::Box const& fill_box,
                     bool const /*overwrite_interior*/,
                     SAMRAI::hier::Transformation const& transformation) const override
    {
        if (!owners_)
            throw std::runtime_error("OwnerFillPattern: no level boxes, setOwners first");

        if constexpr (rank == 0)
        {
            auto const& dst = dynamic_cast<FieldGeometryBase<dimension> const&>(dst_geometry);
            auto const& src = dynamic_cast<FieldGeometryBase<dimension> const&>(src_geometry);
            return std::make_shared<FieldOverlap>(
                destinationBoxes_(dst, src, fill_box, transformation), transformation);
        }
        else
        {
            using Geometry_t = TensorFieldGeometryBase<dimension, rank>;
            auto constexpr N = core::detail::tensor_field_dim_from_rank<rank>();
            auto const& dst  = dynamic_cast<Geometry_t const&>(dst_geometry);
            auto const& src  = dynamic_cast<Geometry_t const&>(src_geometry);
            auto overlaps    = core::for_N_make_array<N>([&](auto i) {
                return std::make_shared<FieldOverlap>(
                    destinationBoxes_(dst[i], src[i], fill_box, transformation), transformation);
            });
            return std::make_shared<TensorFieldOverlap<rank>>(std::move(overlaps));
        }
    }

    std::string const& getPatternName() const override { return s_name_id; }

private:
    // the nodes of dst the copy of src seen through `transformation` owns
    SAMRAI::hier::BoxContainer
    destinationBoxes_(FieldGeometryBase<dimension> const& dst,
                      FieldGeometryBase<dimension> const& src, SAMRAI::hier::Box const& fill_box,
                      SAMRAI::hier::Transformation const& transformation) const
    {
        using dir_t        = SAMRAI::tbox::Dimension::dir_t;
        auto const& shift  = transformation.getOffset();
        auto const& src_id = src.patchBox.getGlobalId();

        SAMRAI::hier::BoxContainer boxes;
        if (src_id == dst.patchBox.getGlobalId()
            and shift == SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension{dimension}))
            return boxes;

        // the nodes of a cell box for this quantity: one more on primal directions
        auto const nodes = [&](SAMRAI::hier::Box box) {
            for (dir_t d = 0; d < box.getDim().getValue(); ++d)
                box.setUpper(d, box.upper(d) + dst.interiorFieldBox().upper(d)
                                    - dst.patchBox.upper(d));
            return box;
        };
        // the nodes a patch may own
        auto const reach = [&](SAMRAI::hier::Box box) {
            if constexpr (band)
                box.grow(SAMRAI::hier::IntVector::getOne(SAMRAI::tbox::Dimension{dimension}));
            return nodes(box);
        };

        SAMRAI::hier::Box src_cells{src.patchBox};
        transformation.transform(src_cells);
        auto const together = dst.ghostFieldBox() * reach(src_cells) * nodes(fill_box);
        if (together.empty())
            return boxes;

        boxes.push_back(together);
        for (auto const& candidate : owners_->candidates)
        {
            auto const theirs = reach(candidate.box);
            if (!theirs.intersects(together))
                continue;
            if constexpr (band) // the nodes of the level are not in the band
                boxes.removeIntersections(nodes(candidate.box));
            if (Owners_t::before(candidate, src_id, shift))
                boxes.removeIntersections(theirs);
        }
        return boxes;
    }

    std::shared_ptr<Owners_t const> owners_;

    // SAMRAI asks the stencil width of any other name (RefineClasses.cpp:204)
    static inline std::string const s_name_id = "BOX_GEOMETRY_FILL_PATTERN";

    SAMRAI::hier::IntVector const& getStencilWidth() override
    {
        throw std::runtime_error("never called");
    }

    std::shared_ptr<SAMRAI::hier::BoxOverlap>
    computeFillBoxesOverlap(SAMRAI::hier::BoxContainer const& /*fill_boxes*/,
                            SAMRAI::hier::BoxContainer const& /*node_fill_boxes*/,
                            SAMRAI::hier::Box const& /*patch_box*/,
                            SAMRAI::hier::Box const& /*data_box*/,
                            SAMRAI::hier::PatchDataFactory const& /*pdf*/) const override
    {
        throw std::runtime_error("OwnerFillPattern: same level schedules only");
    }
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

        auto const& dst_geometry = dynamic_cast<FieldGeometry_t const&>(_dst_geometry);
        auto const& src_geometry = dynamic_cast<FieldGeometry_t const&>(_src_geometry);

        if (isSamePatch(dst_geometry, src_geometry))
            return std::make_shared<FieldOverlap>(SAMRAI::hier::BoxContainer{}, transformation);

        return calculateOverlap(dst_geometry, src_geometry, dst_patch_box, src_mask, fill_box,
                                overwrite_interior, transformation);
    }


    std::shared_ptr<SAMRAI::hier::BoxOverlap> static calculateOverlap(
        auto const& dst_geometry, auto const& src_geometry, SAMRAI::hier::Box const& dst_patch_box,
        SAMRAI::hier::Box const& src_mask, SAMRAI::hier::Box const& fill_box,
        bool const overwrite_interior, SAMRAI::hier::Transformation const& transformation)
    {
        auto const src_ghost_box = [&]() {
            auto const ghost_field_box = phare_box_from<dim>(src_geometry.ghostFieldBox());
            return amr::shift(ghost_field_box, transformation);
        }();

        auto const dst_ghost_box = [&]() {
            auto const ghost_field_box = phare_box_from<dim>(dst_geometry.ghostFieldBox());
            return ghost_field_box;
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

        auto const& dst_geometry = dynamic_cast<TensorFieldGeometry_t const&>(_dst_geometry);
        auto const& src_geometry = dynamic_cast<TensorFieldGeometry_t const&>(_src_geometry);

        if (isSamePatch(dst_geometry, src_geometry))
        {
            auto overlaps = core::for_N_make_array<N>([&](auto /*i*/) {
                return std::make_shared<FieldOverlap>(SAMRAI::hier::BoxContainer{}, transformation);
            });
            return std::make_shared<TensorFieldOverlap<rank_>>(std::move(overlaps));
        }

        auto overlaps = core::for_N_make_array<N>([&](auto i) {
            auto overlap = FieldGhostInterpOverlapFillPattern<Gridlayout_t>::calculateOverlap(
                dst_geometry[i], src_geometry[i], dst_patch_box, src_mask, fill_box,
                overwrite_interior, transformation);

            return std::dynamic_pointer_cast<FieldOverlap>(overlap);
        });
        return std::make_shared<TensorFieldOverlap<rank_>>(std::move(overlaps));
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
