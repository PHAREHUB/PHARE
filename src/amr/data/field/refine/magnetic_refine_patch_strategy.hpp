#ifndef PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
#define PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP

#include "amr/data/field/field_geometry.hpp"
#include "amr/data/field/refine/coarse_cell_round_out.hpp"
#include "amr/data/field/refine/field_refine_patch_strategy.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/data/tensorfield/tensor_field_data_traits.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/utilities/box/amr_box.hpp"

#include "core/utilities/box/box_spans.hpp"
#include "core/data/vector.hpp"
#include "core/data/field/field_box.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/constants.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace PHARE::amr
{

namespace detail
{
    // the untiled tensor field data matching a tiled one, e.g. the phare_scratch_vec_field's
    template<typename T>
    struct untiled_tensor_field_data;

    template<std::size_t rank, typename GridLayoutT, typename Grid_t, typename PhysicalQuantity>
    struct untiled_tensor_field_data<TensorFieldData<rank, GridLayoutT, Grid_t, PhysicalQuantity>>
    {
        using type
            = TensorFieldData<rank, GridLayoutT, typename Grid_t::grid_type, PhysicalQuantity>;
    };

    //! a fine cell of the snapshot region, in the layout's local indices -- never negative
    template<std::size_t dimension>
    using LocalCellKey = std::array<std::uint32_t, dimension>;

    /**
     * @brief The stage-1 divergence snapshot of one touch-up pass: one slot per fine cell of the
     * reconstruction region, NaN meaning "not computed yet".
     *
     * Every correction has to read stage-1 divergences -- recomputing from the live field once a
     * sibling interior face has been written couples the two faces and loses divB exactness -- so
     * each subzone divergence is computed once per pass and reread afterwards. The region is a
     * dense box of fine cells, so the slots are a flat array indexed by the cell's offset in it,
     * and the array is kept between passes: only a region bigger than every earlier one allocates.
     */
    template<typename GridLayoutT>
    class DivScratch
    {
        static constexpr std::size_t dimension = GridLayoutT::dimension;

        using LocalCellKey_t = LocalCellKey<dimension>;

    public:
        //! sizes the snapshot to `region`'s fine cells, in the local indices the corrections
        //! index by, and marks every slot not-computed
        void reset(SAMRAI::hier::Box const& region, GridLayoutT const& layout)
        {
            auto const local = layout.AMRToLocal(phare_box_from<dimension>(region));
            auto const shape = local.shape();

            for (std::size_t d = 0; d < dimension; ++d)
            {
                lower_[d] = static_cast<std::uint32_t>(local.lower[d]);
                shape_[d] = static_cast<std::uint32_t>(shape[d]);
            }

            values_.assign(local.size(), unset_);
        }

        NO_DISCARD double& operator()(core::Point<int, dimension> const& cell)
        {
            LocalCellKey_t key;
            for (std::size_t d = 0; d < dimension; ++d)
                key[d] = static_cast<std::uint32_t>(cell[d]);

            return values_[offset_(key)];
        }

    private:
        static constexpr double unset_ = std::numeric_limits<double>::quiet_NaN();

        // row-major offset of a fine cell in the region. The corrections only ever reach cells of
        // the coarse cell they correct an interior face of, so a region of whole coarse cells
        // (reconstructionRegion) contains every cell they index -- the assert catches a caller
        // that reset the snapshot to something else.
        NO_DISCARD std::size_t offset_(LocalCellKey_t const& cell) const
        {
            std::size_t offset = 0;

            for (std::size_t d = 0; d < dimension; ++d)
            {
                // unsigned: a cell under lower_ wraps around and trips the same bound
                assert(cell[d] - lower_[d] < shape_[d]);

                offset = offset * std::size_t{shape_[d]} + std::size_t{cell[d] - lower_[d]};
            }

            return offset;
        }

        core::Point<std::uint32_t, dimension> lower_;
        core::Point<std::uint32_t, dimension> shape_;
        std::vector<double> values_;
    };
} // namespace detail


/**
 * @brief Stage 2 of the Balsara divergence-free magnetic prolongation: the cross-component
 *        divB touch-up. Order-independent.
 *
 * Stage 1 (the fill-all composite kernel, magnetic_composite_refiner.hpp) fills every fine face of
 * each B component from that component's own coarse faces, and makes no divB claim. This stage
 * adds a closed-form min-norm correction to the 2d interior faces of each coarse zone, equalizing
 * the divergences of the 2^d fine cells the zone splits into -- its subzones. They all become the
 * zone's transported divergence q0, so a discretely div-free coarse field (q0 = 0) prolongs
 * div-free to roundoff whatever the stage-1 order.
 *
 * Balsara, D. S., Samantaray, S. & Subramanian, S. (2023), "Efficient WENO-Based Prolongation
 * Strategies for Divergence-Preserving Vector Fields", Commun. Appl. Math. Comput.,
 * doi:10.1007/s42967-021-00182-x. The correction of an interior face is a weighted sum of the
 * zone's subzone divergence deficits — own subzone weighted highest, each farther
 * (Hamming-distance) neighbour less. The per-dimension weights are the flux-variable ones,
 * 1/2 (1D), [3,1]/8 (2D) and [7,2,2,1]/24 (3D): the density weights ([·]/16, [·]/48) rescaled by
 * the refinement ratio 2.
 *
 * The correction is added to the stage-1 face value rather than replacing it, so a higher-order
 * stage-1 interior fill survives the touch-up.
 */
template<IsVecFieldData VecFieldDataT, typename BoundaryManagerT>
class MagneticRefinePatchStrategy : public FieldRefinePatchStrategy<VecFieldDataT, BoundaryManagerT>
{
public:
    using geometry_type       = VecFieldDataT::Geometry;
    using field_geometry_type = geometry_type::FieldGeometry_t;
    using gridlayout_type     = VecFieldDataT::gridlayout_type;
    using Super               = FieldRefinePatchStrategy<VecFieldDataT, BoundaryManagerT>;

    static constexpr std::size_t dimension = VecFieldDataT::dimension;
    static constexpr std::size_t N         = VecFieldDataT::N;


    using Super::data_id_;

    using DivScratch = detail::DivScratch<gridlayout_type>;
    using Point_t    = core::Point<int, dimension>;

    /**
     * @brief Construct the magnetic refinement strategy.
     * @param boundary_manager Manager handling boundary conditions.
     */
    MagneticRefinePatchStrategy(BoundaryManagerT& boundaryManager)
        : Super(boundaryManager)
    {
    }

    // the untiled phare_scratch_vec_field, tiled B gathers its touch-up region into
    void registerScratchID(int const scratch_id) { scratch_id_ = scratch_id; }

    void preprocessRefine(SAMRAI::hier::Patch&, SAMRAI::hier::Patch const&,
                          SAMRAI::hier::Box const&, SAMRAI::hier::IntVector const&) override
    {
    }

    // Always 0: the touch-up reads only fine faces the stage-1 gather already wrote, no coarse
    // data. SAMRAI provisions max(this, every registered refine operator), so 0 narrows nothing.
    SAMRAI::hier::IntVector
    getRefineOpStencilWidth(SAMRAI::tbox::Dimension const& dim) const override
    {
        return SAMRAI::hier::IntVector::getZero(dim);
    }

    /**
     * @brief The region the touch-up over fill box `fine_box` must run over: `fine_box` rounded
     * out to whole coarse cells, clipped to the allocation.
     *
     * The touch-up reaches exactly one coarse cell, so it is self-consistent only over a union of
     * whole coarse cells (coarse_cell_round_out.hpp). SAMRAI fill boxes are not: a recursive
     * schedule's coarse-interpolation temporary plus its one-cell ring always cuts a coarse cell
     * in half, leaving the far shared faces unwritten. Stage 1 fills them from the coarse level,
     * except where a copy of same-level data already did (copiesWholeCoarseCells).
     *
     * Round-out grows each edge by one index, so the region stays inside the allocation iff
     * field_ghost_width >= fill_ring + 1. That holds everywhere we support, but with zero slack:
     * the ring is 1 and field_ghost_width (ghost_width_calculator.hpp) bottoms out at exactly 2.
     * If the clip ever bit it would halve a coarse cell, and a half-covered cell has no well-posed
     * reconstruction — some of its inputs are faces nothing ever wrote — hence the throw.
     */
    static SAMRAI::hier::Box reconstructionRegion(SAMRAI::hier::Box const& fine_box,
                                                  SAMRAI::hier::Box const& ghost_box)
    {
        auto const region = roundCellBoxOutToCoarseCells<dimension>(fine_box) * ghost_box;

        if (!isWholeCoarseCells<dimension>(region))
            throw std::runtime_error(
                "magnetic prolongation region is not a union of whole coarse cells: fill box "
                + to_str(fine_box) + ", rounded out and clipped to the allocated "
                + to_str(ghost_box) + ", gives " + to_str(region)
                + ", which half-covers a coarse cell. The field ghost width must be at least the "
                  "fill ring plus one.");

        return region;
    }




    // Run over the whole-coarse-cell-rounded region so every subzoneDiv_ read lands on a fine
    // face the stage-1 refinement actually wrote.
    void postprocessRefine(SAMRAI::hier::Patch& fine, SAMRAI::hier::Patch const& coarse,
                           SAMRAI::hier::Box const& fine_box,
                           SAMRAI::hier::IntVector const& ratio) override
    {
        Super::assertIDsSet();

        auto& fields = VecFieldDataT::getFields(fine, data_id_);

        auto const layout = PHARE::amr::layoutFromPatch<gridlayout_type>(fine);
        auto const region
            = reconstructionRegion(fine_box, fine.getPatchData(data_id_)->getGhostBox());

        using Field_t = std::decay_t<decltype(fields[0])>;
        if constexpr (core::is_field_tile_set_v<Field_t>)
            touchUpTiles(fine, fields, layout, region);
        else
            touchUpInteriorFaces(fields, layout, region, scratch_);
    }


    // the touch-up runs on patch shaped views of the untiled phare_scratch_vec_field (or of tmp on
    // the SAMRAI temporary patches which have none), holding the owning tiles' values over the
    // region field boxes, which are then copied into every tile's ghost copies
    void touchUpTiles(SAMRAI::hier::Patch const& fine, auto& fields, gridlayout_type const& layout,
                      SAMRAI::hier::Box const& region)
    {
        using View_t       = std::decay_t<decltype(fields[0])>::field_type;
        using ScratchDataT = detail::untiled_tensor_field_data<VecFieldDataT>::type;

        if (scratch_id_ < 0)
            throw std::runtime_error(
                "MagneticRefinePatchStrategy: registerScratchID was not called before use");

        auto const regionLayout = geometry_type::layoutFromBox(region, layout);
        auto const shapes       = core::for_N_make_array<N>(
            [&](auto i) { return layout.allocSize(fields[i].physicalQuantity()); });

        auto const data = [&]() {
            std::array<value_type*, N> ptrs;
            if (fine.checkAllocated(scratch_id_))
            {
                auto& scratch = ScratchDataT::getFields(fine, scratch_id_);
                for (std::size_t i = 0; i < N; ++i)
                {
                    if (core::product(shapes[i]) > scratch[i].size())
                        throw std::runtime_error(
                            "MagneticRefinePatchStrategy: scratch field too small");
                    ptrs[i] = scratch[i].data();
                }
                return ptrs;
            }
            auto& buffer = tmp.get_no_copy(
                core::sum_from(shapes, [](auto const& shape) { return core::product(shape); }));
            for (std::size_t i = 0, offset = 0; i < N; offset += core::product(shapes[i++]))
                ptrs[i] = buffer.data() + offset;
            return ptrs;
        }();

        auto const boxes = core::for_N_make_array<N>([&](auto i) {
            auto const qty              = fields[i].physicalQuantity();
            using PhysicalQuantity      = std::decay_t<decltype(qty)>;
            auto const region_field_box = phare_box_from<dimension>(
                FieldGeometry<gridlayout_type, PhysicalQuantity>::toFieldBox(region, qty,
                                                                             regionLayout));
            return *(grow(region_field_box, 1) * layout.AMRGhostBoxFor(qty));
        });

        auto views = core::for_N_make_array<N>(
            [&](auto i) { return View_t{fields[i].physicalQuantity(), data[i], shapes[i]}; });

        for (std::size_t i = 0; i < N; ++i)
            core::gather_owned(fields[i], boxes[i], views[i],
                               layout.AMRGhostBoxFor(fields[i].physicalQuantity()).lower);

        touchUpInteriorFaces(views, layout, region, scratch_);

        for (std::size_t i = 0; i < N; ++i)
            core::scatter_to_tiles(views[i],
                                   layout.AMRGhostBoxFor(fields[i].physicalQuantity()).lower,
                                   boxes[i], fields[i]);
    }


    /**
     * @brief Add the divergence-equalizing correction over `region`, which must already be
     * whole-coarse-cell rounded (reconstructionRegion): every read here lands on a fine face of
     * the coarse cell being corrected, and those exist only if the cell is wholly inside.
     */
    static void touchUpInteriorFaces(auto& fields, gridlayout_type const& layout,
                                     SAMRAI::hier::Box const& region, DivScratch& scratch)
    {
        auto const regionLayout = geometry_type::layoutFromBox(region, layout);

        // One stage-1 snapshot per pass (DivScratch): every correction must read stage-1
        // divergences, and the corrections write the faces those divergences are computed from.
        scratch.reset(region, layout);

        core::for_N<dimension>([&](auto ic) {
            constexpr auto normal = ic();

            auto const faces = phare_box_from<dimension>(field_geometry_type::toFieldBox(
                region, fields[normal].physicalQuantity(), regionLayout));

            core::for_each_point(faces, [&](auto const& face) {
                correctFace<normal>(scratch, fields, layout, face);
            });
        });
    }


    /**
     * @brief Add the divergence-equalizing correction to the face of component `normal` at AMR
     * index `amrFace`, if that face is interior to its coarse cell.
     *
     * Faces at an even index along the normal coincide with coarse faces and keep their stage-1
     * value. An odd face separates fine cells `face - e_normal` and `face`. The coarse cell holds
     * 2^(d-1) such cell pairs across its normal midplane: the face's own, one step toward the
     * sibling along each transverse direction, and in 3D the diagonal one. The correction is the
     * weighted sum of their divergence jumps, own pair weighted highest ([1]/2, [3,1]/8,
     * [7,2,2,1]/24), times the mesh size along the normal.
     */
    template<std::size_t normal>
    static void correctFace(DivScratch& scratch, auto& fields, gridlayout_type const& layout,
                            Point_t const& amrFace)
    {
        if (!isOddIndex(amrFace[normal]))
            return;

        auto const face = layout.AMRToLocal(amrFace).as_signed();
        auto const& D   = layout.meshSize();
        auto& B         = fields[normal];

        // divergence jump across the face, in the cell pair shifted transversally by `shift`
        auto jump = [&](Point_t const& shift) {
            auto const high = face + shift;
            return subzoneDiv_(scratch, fields, D, high)
                   - subzoneDiv_(scratch, fields, D, high - unit_(normal));
        };

        constexpr auto transverse = transverseDirections_<normal>();
        auto const own            = Point_t{};

        if constexpr (dimension == 1)
        {
            B(face) += D[normal] * jump(own) / 2.0;
        }
        else if constexpr (dimension == 2)
        {
            auto const t = towardSibling_(amrFace, transverse[0]);

            B(face) += D[normal] * (3.0 * jump(own) + jump(t)) / 8.0;
        }
        else
        {
            auto const t1 = towardSibling_(amrFace, transverse[0]);
            auto const t2 = towardSibling_(amrFace, transverse[1]);

            B(face) += D[normal]
                       * (7.0 * jump(own) + 2.0 * jump(t1) + 2.0 * jump(t2) + jump(t1 + t2)) / 24.0;
        }
    }


private:
    // the directions other than `normal`, in increasing order
    template<std::size_t normal>
    static constexpr auto transverseDirections_()
    {
        std::array<std::size_t, dimension - 1> directions{};
        for (std::size_t d = 0, k = 0; d < dimension; ++d)
            if (d != normal)
                directions[k++] = d;
        return directions;
    }

    static Point_t unit_(std::size_t const direction)
    {
        Point_t e{};
        e[direction] = 1;
        return e;
    }

    // one fine cell along `direction` toward the sibling inside the same coarse cell
    static Point_t towardSibling_(Point_t const& amrFace, std::size_t const direction)
    {
        Point_t shift{};
        shift[direction] = isOddIndex(amrFace[direction]) ? -1 : +1;
        return shift;
    }

    // 1/D-weighted divergence of the fine cell at local index `cell`: the sum over directions of
    // (high face − low face)/D_d. Held in the pass snapshot so a correction still reads the
    // stage-1 value after sibling faces are written.
    static double subzoneDiv_(DivScratch& scratch, auto const& fields, auto const& D,
                              Point_t const& cell)
    {
        auto& div = scratch(cell);

        if (std::isnan(div))
        {
            auto faceDifference = [&](std::size_t const d) {
                return (fields[d](cell + unit_(d)) - fields[d](cell)) / D[d];
            };

            div = faceDifference(0);
            for (std::size_t d = 1; d < dimension; ++d)
                div += faceDifference(d);
        }

        return div;
    }

    int scratch_id_ = -1;
    DivScratch scratch_;

    using value_type = VecFieldDataT::value_type;
    static inline core::MinimizingVector<value_type> tmp; // SAMRAI temporary patches
};

} // namespace PHARE::amr

#endif // PHARE_AMR_MAGNETIC_REFINE_PATCH_STRATEGY_HPP
