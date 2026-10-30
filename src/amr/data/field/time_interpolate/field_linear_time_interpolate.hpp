#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "phare_mpi.hpp" // IWYU pragma: keep


#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/field/field_overlap.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "amr/data/tensorfield/tensor_field_overlap.hpp"

#include <SAMRAI/hier/TimeInterpolateOperator.h>

#include <tuple>



namespace PHARE::amr
{


template<typename Dst>
void linear_time_interpolate(Dst& fieldDest, auto const& fieldSrcOld, auto const& fieldSrcNew,
                             auto&&... args)
{
    auto const& [localDestBox, localSrcBox, alpha] = std::forward_as_tuple(args...);
    auto const lclDstBox                           = phare_box_from<Dst::dimension>(localDestBox);
    auto const lclSrcBox                           = phare_box_from<Dst::dimension>(localSrcBox);

    auto src_it = lclSrcBox.begin();
    auto dst_it = lclDstBox.begin();

    for (; dst_it != lclDstBox.end(); ++src_it, ++dst_it)
        fieldDest(*dst_it) = (1. - alpha) * fieldSrcOld(*src_it) + alpha * fieldSrcNew(*src_it);
}


/**
 * @brief Call `interpolate(localDestBox, localSrcBox)` over every destination box of `overlap`.
 *
 * Not over SAMRAI's `where` (the source mask): a rounded B overlap
 * (FieldGeometry::calculateOverlap) is larger. Under an offset the destination is a temporary in
 * source space, hence the inverse transform, as in FieldData::packStream.
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


template<typename GridLayoutT, typename FieldT>
class FieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
    static std::size_t constexpr dim = GridLayoutT::dimension;
    static_assert(dim > 0 && dim <= 3);

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutT, FieldT>;

public:
    using GridLayoutImpl = GridLayoutT::implT;

    FieldLinearTimeInterpolate()
        : SAMRAI::hier::TimeInterpolateOperator{"FieldLinearTimeInterpolate"}
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

        auto const& interpTime = fieldDataDest.getTime();
        auto const& oldTime    = fieldDataSrcOld.getTime();
        auto const& newTime    = fieldDataSrcNew.getTime();
        auto const& alpha      = (interpTime - oldTime) / (newTime - oldTime);

        auto const& fieldSrcOld = fieldDataSrcOld.field;
        auto const& fieldSrcNew = fieldDataSrcNew.field;
        auto& fieldDest         = fieldDataDest.field;

        auto const& layout = fieldDataDest.gridLayout;
        auto qty           = fieldDest.physicalQuantity();

        auto const ghostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            fieldDataDest.getGhostBox(), qty, layout);

        auto srcGhostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            fieldDataSrcNew.getGhostBox(), qty, fieldDataSrcNew.gridLayout);

        interpolateOverOverlap(dynamic_cast<FieldOverlap const&>(overlap), ghostBox, srcGhostBox,
                               [&](auto const& localDestBox, auto const& localSrcBox) {
                                   linear_time_interpolate(fieldDest, fieldSrcOld, fieldSrcNew,
                                                           localDestBox, localSrcBox, alpha);
                               });
    }
};


template<std::size_t rank, typename GridLayoutT, typename FieldT, typename PhysicalQuantity>
class TensorFieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
    static std::size_t constexpr dim = GridLayoutT::dimension;
    static_assert(dim > 0 && dim <= 3);

    using TensorFieldDataT         = TensorFieldData<rank, GridLayoutT, FieldT, PhysicalQuantity>;
    static constexpr std::size_t N = TensorFieldDataT::N;

public:
    using GridLayoutImpl = typename GridLayoutT::implT;

    TensorFieldLinearTimeInterpolate()
        : SAMRAI::hier::TimeInterpolateOperator{"TensorFieldLinearTimeInterpolate"}
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

            interpolateOverOverlap(*tensorOverlap[c], ghostBox, srcGhostBox,
                                   [&](auto const& localDestBox, auto const& localSrcBox) {
                                       linear_time_interpolate(fieldDests[c], fieldSrcOlds[c],
                                                               fieldSrcNews[c], localDestBox,
                                                               localSrcBox, alpha);
                                   });
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename PhysicalQuantity>
using VecFieldLinearTimeInterpolate
    = TensorFieldLinearTimeInterpolate</*rank=*/1, GridLayoutT, FieldT, PhysicalQuantity>;


} // namespace PHARE::amr

#endif
