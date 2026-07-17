#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/grid/grid_tiles.hpp"

#include "amr/data/field/field_data.hpp"
#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"

#include <SAMRAI/hier/TimeInterpolateOperator.h>

#include <tuple>



namespace PHARE::amr
{


template<typename Dst> // ONLY TAKES PHARE::core::Box!
void linear_time_interpolate(Dst& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                             auto&&... args)
{
    auto const& [lclDstBox, lclSrcBox, alpha] = std::forward_as_tuple(args...);

    auto src_it = lclSrcBox.begin();
    auto dst_it = lclDstBox.begin();

    for (; dst_it != lclDstBox.end(); ++src_it, ++dst_it)
        fieldDest(*dst_it) = (1. - alpha) * fieldSrcOld(*src_it) + alpha * fieldSrcNew(*src_it);
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

    void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& where,
                         SAMRAI::hier::BoxOverlap const& /*overlap*/,
                         SAMRAI::hier::PatchData const& srcDataOld,
                         SAMRAI::hier::PatchData const& srcDataNew) const override
    {
        auto& fieldDataDest = dynamic_cast<FieldDataT&>(destData);

        auto const& fieldDataSrcOld = dynamic_cast<FieldDataT const&>(srcDataOld);
        auto const& fieldDataSrcNew = dynamic_cast<FieldDataT const&>(srcDataNew);

        auto const& fieldSrcOld = fieldDataSrcOld.field;
        auto const& fieldSrcNew = fieldDataSrcNew.field;
        auto& fieldDest         = fieldDataDest.field;

        auto const& layout     = fieldDataDest.gridLayout;
        auto const whereLayout = FieldGeometry_t::layoutFromBox(where, layout);

        auto const interpTime = fieldDataDest.getTime();
        auto const oldTime    = fieldDataSrcOld.getTime();
        auto const newTime    = fieldDataSrcNew.getTime();
        auto const alpha      = (interpTime - oldTime) / (newTime - oldTime);
        auto const qty        = fieldDest.physicalQuantity();

        auto const dstGhostBox = phare_box_from<dim>(
            FieldGeometry_t::toFieldBox(fieldDataDest.getGhostBox(), qty, layout));

        auto const interpolateBox
            = phare_box_from<dim>(FieldGeometry_t::toFieldBox(where, qty, whereLayout));

        auto const srcGhostBox = phare_box_from<dim>(FieldGeometry_t::toFieldBox(
            fieldDataSrcNew.getGhostBox(), qty, fieldDataSrcNew.gridLayout));

        interpolate(fieldDest, dstGhostBox, fieldSrcOld, fieldSrcNew, srcGhostBox, interpolateBox,
                    qty, alpha);
    }

    void timeInterpolateField(auto& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                              auto const qty, auto const alpha, auto const ghostBox,
                              auto const srcGhostBox, auto const finalBox) const
    {
        auto const localDestBox = AMRToLocal(finalBox, ghostBox);
        auto const localSrcBox  = AMRToLocal(finalBox, srcGhostBox);

        linear_time_interpolate( //
            fieldDest, fieldSrcOld, fieldSrcNew, localDestBox, localSrcBox, alpha);
    }



    void interpolate(auto& fieldDest, auto& dstGhostBox, auto& fieldSrcOld, auto& fieldSrcNew,
                     auto& srcGhostBox, auto& box, auto& qty, auto& alpha) const
    {
        if constexpr (core::is_field_tile_set_v<FieldT>)
        {
            for (std::size_t tidx = 0; tidx < fieldSrcOld().size(); ++tidx)
            {
                auto& srcTileOld   = fieldSrcOld()[tidx];
                auto& srcTileNew   = fieldSrcNew()[tidx];
                auto const src_box = srcTileOld.ghost_box();

                if (auto const src_overlap = src_box * srcGhostBox)
                    for (auto& dst_tile : fieldDest())
                    {
                        auto const dst_box = dst_tile.ghost_box();

                        if (auto const dst_overlap = dst_box * dstGhostBox)
                        {
                            auto const finalBox = box * (*(*dst_overlap * *src_overlap));

                            timeInterpolateField(dst_tile(), srcTileOld(), srcTileNew(), qty, alpha,
                                                 dst_box, src_box, *finalBox);
                        }
                    }
            }
        }
        else
        {
            auto const finalBox = box * dstGhostBox;
            timeInterpolateField(fieldDest, fieldSrcOld, fieldSrcNew, qty, alpha, dstGhostBox,
                                 srcGhostBox, *finalBox);
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
    using GridLayoutImpl = typename GridLayoutT::implT;

    TensorFieldLinearTimeInterpolate()
        : Super{"TensorFieldLinearTimeInterpolate"}
    {
    }


    virtual ~TensorFieldLinearTimeInterpolate() = default;


    void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& where,
                         SAMRAI::hier::BoxOverlap const& /*overlap*/,
                         SAMRAI::hier::PatchData const& srcDataOld,
                         SAMRAI::hier::PatchData const& srcDataNew) const override
    {
        auto& fieldDataDest = dynamic_cast<TensorFieldDataT&>(destData);

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

            auto const& whereLayout = FieldGeometry_t::layoutFromBox(where, layout);
            auto const& interpolateBox
                = phare_box_from<dim>(FieldGeometry_t::toFieldBox(where, qty, whereLayout));
            auto const& dstGhostBox = phare_box_from<dim>(
                FieldGeometry_t::toFieldBox(fieldDataDest.getGhostBox(), qty, layout));

            auto const& srcGhostBox = phare_box_from<dim>(FieldGeometry_t::toFieldBox(
                fieldDataSrcNew.getGhostBox(), qty, fieldDataSrcNew.gridLayout));

            Super::interpolate(fieldDests[c], dstGhostBox, fieldSrcOlds[c], fieldSrcNews[c],
                               srcGhostBox, interpolateBox, qty, alpha);
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename PhysicalQuantity>
using VecFieldLinearTimeInterpolate
    = TensorFieldLinearTimeInterpolate</*rank=*/1, GridLayoutT, FieldT, PhysicalQuantity>;


} // namespace PHARE::amr

#endif
