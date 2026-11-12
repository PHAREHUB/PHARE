#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/grid/grid_tiles.hpp"

#include "amr/data/field/field_data.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "core/utilities/box/box_spans.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/field/field_overlap.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"

#include <SAMRAI/hier/TimeInterpolateOperator.h>

#include <tuple>



namespace PHARE::amr
{


template<typename Dst> // ONLY TAKES PHARE::core::Box!
void linear_time_interpolate(Dst& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                             auto&&... args)
{
    auto const& [lclDstBox, lclSrcBox, alpha] = std::forward_as_tuple(args...);

    core::for_each_point_zip(lclDstBox, lclSrcBox, [&](auto const& dst, auto const& src) {
        fieldDest(dst) = (1. - alpha) * fieldSrcOld(src) + alpha * fieldSrcNew(src);
    });
}


/**
 * @brief Call `interpolate(localDestBox, localSrcBox)` over every destination box of `overlap`.
 *
 * The overlap is used rather than SAMRAI's `where` (the source mask), because a B overlap rounded
 * out to whole coarse cells (FieldGeometry::calculateOverlap) extends past it. Under a periodic
 * offset the destination data is a temporary located in source space, so the destination boxes
 * are inverse-transformed, as in FieldData::packStream.
 */
void interpolateOverOverlap(FieldOverlap const& overlap, SAMRAI::hier::Box const& ghostBox,
                            SAMRAI::hier::Box const& srcGhostBox, auto&& interpolate)
{
    for (auto box : overlap.getDestinationBoxContainer())
    {
        overlap.getTransformation().inverseTransform(box);
        auto const finalBox = box * ghostBox * srcGhostBox;
        if (!finalBox.empty())
            interpolate(AMRToLocal(finalBox, ghostBox), AMRToLocal(finalBox, srcGhostBox));
    }
}


/**
 * @brief Call `fn(box)` with every destination box of `overlap`, inverse-transformed as in
 *        interpolateOverOverlap, for tiled fields to intersect with their tiles.
 */
void forOverlapBoxes(FieldOverlap const& overlap, auto&& fn)
{
    for (auto box : overlap.getDestinationBoxContainer())
    {
        overlap.getTransformation().inverseTransform(box);
        fn(box);
    }
}


template<typename GridLayoutT, typename FieldT>
class FieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
    static std::size_t constexpr dim = GridLayoutT::dimension;
    static_assert(dim > 0 && dim <= 3);

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutT, FieldT>;
    using FieldGeometry_t  = FieldGeometry<GridLayoutT, PhysicalQuantity>;
    bool static constexpr withGhost{true};

    auto static toFieldBox(auto&&... args) { return FieldGeometry_t::toFieldBox(args...); }

protected:
    FieldLinearTimeInterpolate(std::string const& name)
        : SAMRAI::hier::TimeInterpolateOperator{name}
    {
    }

public:
    using GridLayoutImpl = GridLayoutT::implT;

    FieldLinearTimeInterpolate()
        : FieldLinearTimeInterpolate{"FieldLinearTimeInterpolate"}
    {
    }


    virtual ~FieldLinearTimeInterpolate() = default;

    void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& /*where*/,
                         SAMRAI::hier::BoxOverlap const& overlap,
                         SAMRAI::hier::PatchData const& srcDataOld,
                         SAMRAI::hier::PatchData const& srcDataNew) const override
    {
        auto& fieldDataDest = dynamic_cast<FieldDataT&>(destData);

        auto const& fieldDataSrcOld = dynamic_cast<FieldDataT const&>(srcDataOld);
        auto const& fieldDataSrcNew = dynamic_cast<FieldDataT const&>(srcDataNew);

        auto const& fieldSrcOld = fieldDataSrcOld.field;
        auto const& fieldSrcNew = fieldDataSrcNew.field;
        auto& fieldDest         = fieldDataDest.field;

        auto const& layout    = fieldDataDest.gridLayout;
        auto const interpTime = fieldDataDest.getTime();
        auto const oldTime    = fieldDataSrcOld.getTime();
        auto const newTime    = fieldDataSrcNew.getTime();
        auto const alpha      = (interpTime - oldTime) / (newTime - oldTime);
        auto const qty        = fieldDest.physicalQuantity();

        auto const ghostBox = toFieldBox(fieldDataDest.getGhostBox(), qty, layout);
        auto const srcGhostBox
            = toFieldBox(fieldDataSrcNew.getGhostBox(), qty, fieldDataSrcNew.gridLayout);

        interpolate(fieldDest, fieldSrcOld, fieldSrcNew, dynamic_cast<FieldOverlap const&>(overlap),
                    ghostBox, srcGhostBox, alpha);
    }

    void timeInterpolateField(auto& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                              auto const alpha, auto const ghostBox, auto const srcGhostBox,
                              auto const finalBox) const
    {
        auto const localDestBox = AMRToLocal(finalBox, ghostBox);
        auto const localSrcBox  = AMRToLocal(finalBox, srcGhostBox);

        linear_time_interpolate( //
            fieldDest, fieldSrcOld, fieldSrcNew, localDestBox, localSrcBox, alpha);
    }


    void interpolate(auto& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                     FieldOverlap const& overlap, SAMRAI::hier::Box const& ghostBox,
                     SAMRAI::hier::Box const& srcGhostBox, double const alpha) const
    {
        if constexpr (core::is_field_tile_set_v<FieldT>)
        {
            auto const dstGhostBox = phare_box_from<dim>(ghostBox);
            auto const srcGhost    = phare_box_from<dim>(srcGhostBox);

            forOverlapBoxes(overlap, [&](auto const& amr_box) {
                auto const box = phare_box_from<dim>(amr_box);

                for (std::size_t tidx = 0; tidx < fieldSrcOld().size(); ++tidx)
                {
                    auto& srcTileOld   = fieldSrcOld()[tidx];
                    auto& srcTileNew   = fieldSrcNew()[tidx];
                    auto const src_box = srcTileOld.ghost_box();

                    if (auto const src_overlap = src_box * srcGhost)
                        for (auto& dst_tile : fieldDest())
                        {
                            auto const dst_box = dst_tile.ghost_box();

                            if (auto const dst_overlap = dst_box * dstGhostBox)
                                if (auto const tile_overlap = *dst_overlap * *src_overlap)
                                    if (auto const finalBox = box * *tile_overlap)
                                        timeInterpolateField(dst_tile(), srcTileOld(),
                                                             srcTileNew(), alpha, dst_box, src_box,
                                                             *finalBox);
                        }
                }
            });
        }
        else
        {
            interpolateOverOverlap(overlap, ghostBox, srcGhostBox,
                                   [&](auto const& localDestBox, auto const& localSrcBox) {
                                       // SAMRAI boxes, linear_time_interpolate takes core::Box
                                       linear_time_interpolate(
                                           fieldDest, fieldSrcOld, fieldSrcNew,
                                           phare_box_from<dim>(localDestBox),
                                           phare_box_from<dim>(localSrcBox), alpha);
                                   });
        }
    }
};


template<std::size_t rank, typename GridLayoutT, typename FieldT, typename PhysicalQuantity>
class TensorFieldLinearTimeInterpolate : public FieldLinearTimeInterpolate<GridLayoutT, FieldT>
{
    static std::size_t constexpr dim = GridLayoutT::dimension;
    static_assert(dim > 0 && dim <= 3);

    using Super                    = FieldLinearTimeInterpolate<GridLayoutT, FieldT>;
    using TensorFieldDataT         = TensorFieldData<rank, GridLayoutT, FieldT, PhysicalQuantity>;
    static constexpr std::size_t N = TensorFieldDataT::N;

public:
    using GridLayoutImpl = GridLayoutT::implT;

    TensorFieldLinearTimeInterpolate()
        : Super{"TensorFieldLinearTimeInterpolate"}
    {
    }


    virtual ~TensorFieldLinearTimeInterpolate() = default;


    void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& /*where*/,
                         SAMRAI::hier::BoxOverlap const& overlap,
                         SAMRAI::hier::PatchData const& srcDataOld,
                         SAMRAI::hier::PatchData const& srcDataNew) const override
    {
        auto& fieldDataDest       = dynamic_cast<TensorFieldDataT&>(destData);
        auto const& tensorOverlap = dynamic_cast<TensorFieldOverlap<rank> const&>(overlap);

        auto const& fieldDataSrcOld = dynamic_cast<TensorFieldDataT const&>(srcDataOld);
        auto const& fieldDataSrcNew = dynamic_cast<TensorFieldDataT const&>(srcDataNew);

        auto const& interpTime   = fieldDataDest.getTime();
        auto const& oldTime      = fieldDataSrcOld.getTime();
        auto const& newTime      = fieldDataSrcNew.getTime();
        auto const& alpha        = (interpTime - oldTime) / (newTime - oldTime);
        auto const& fieldSrcOlds = fieldDataSrcOld.grids;
        auto const& fieldSrcNews = fieldDataSrcNew.grids;
        auto& fieldDests         = fieldDataDest.grids;
        auto const& layout       = fieldDataDest.gridLayout;

        for (std::uint16_t c = 0; c < N; ++c)
        {
            auto const& qty       = fieldDests[c].physicalQuantity();
            using FieldGeometry_t = FieldGeometry<GridLayoutT, std::decay_t<decltype(qty)>>;

            auto const& ghostBox
                = FieldGeometry_t::toFieldBox(fieldDataDest.getGhostBox(), qty, layout);
            auto const& srcGhostBox = FieldGeometry_t::toFieldBox(fieldDataSrcNew.getGhostBox(),
                                                                  qty, fieldDataSrcNew.gridLayout);

            Super::interpolate(fieldDests[c], fieldSrcOlds[c], fieldSrcNews[c], *tensorOverlap[c],
                               ghostBox, srcGhostBox, alpha);
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename PhysicalQuantity>
using VecFieldLinearTimeInterpolate
    = TensorFieldLinearTimeInterpolate</*rank=*/1, GridLayoutT, FieldT, PhysicalQuantity>;


} // namespace PHARE::amr

#endif
