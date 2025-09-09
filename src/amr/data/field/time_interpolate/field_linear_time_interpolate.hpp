#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/tensorfield/tensor_field_data.hpp"

#include <SAMRAI/hier/TimeInterpolateOperator.h>

#include <tuple>



namespace PHARE::amr
{


template<typename Dst>
void linear_time_interpolate(Dst& fieldDest, auto& fieldSrcOld, auto& fieldSrcNew, auto&&... args)
{
    auto const& [localDestBox, localSrcBox, alpha] = std::forward_as_tuple(args...);
    auto const lclDstBox                           = phare_box_from<Dst::dimension>(localDestBox);
    auto const lclSrcBox                           = phare_box_from<Dst::dimension>(localSrcBox);

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

public:
    using GridLayoutImpl = typename GridLayoutT::implT;

    FieldLinearTimeInterpolate()
        : SAMRAI::hier::TimeInterpolateOperator{"FieldLinearTimeInterpolate"}
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

        auto const& interpTime = fieldDataDest.getTime();
        auto const& oldTime    = fieldDataSrcOld.getTime();
        auto const& newTime    = fieldDataSrcNew.getTime();
        auto const& alpha      = (interpTime - oldTime) / (newTime - oldTime);

        auto const& fieldSrcOld = fieldDataSrcOld.field;
        auto const& fieldSrcNew = fieldDataSrcNew.field;
        auto& fieldDest         = fieldDataDest.field;

        auto const& layout = fieldDataDest.gridLayout;
        auto const whereLayout
            = FieldGeometry<GridLayoutT, PhysicalQuantity>::layoutFromBox(where, layout);

        auto qty = fieldDest.physicalQuantity();
        auto const interpolateBox
            = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(where, qty, whereLayout);

        auto const ghostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            fieldDataDest.getGhostBox(), qty, layout);

        auto const finalBox = interpolateBox * ghostBox;

        auto srcGhostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            fieldDataSrcNew.getGhostBox(), qty, fieldDataSrcNew.gridLayout);

        auto const localDestBox = AMRToLocal(finalBox, ghostBox);
        auto const localSrcBox  = AMRToLocal(finalBox, srcGhostBox);

        linear_time_interpolate( //
            fieldDest, fieldSrcOld, fieldSrcNew, localDestBox, localSrcBox, alpha);
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
        : SAMRAI::hier::TimeInterpolateOperator{"FieldLinearTimeInterpolate"}
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

            auto const& whereLayout    = FieldGeometry_t::layoutFromBox(where, layout);
            auto const& interpolateBox = FieldGeometry_t::toFieldBox(where, qty, whereLayout);
            auto const& ghostBox
                = FieldGeometry_t::toFieldBox(fieldDataDest.getGhostBox(), qty, layout);
            auto const& finalBox     = interpolateBox * ghostBox;
            auto const& srcGhostBox  = FieldGeometry_t::toFieldBox(fieldDataSrcNew.getGhostBox(),
                                                                   qty, fieldDataSrcNew.gridLayout);
            auto const& localDestBox = AMRToLocal(finalBox, ghostBox);
            auto const& localSrcBox  = AMRToLocal(finalBox, srcGhostBox);

            linear_time_interpolate( //
                fieldDests[c], fieldSrcOlds[c], fieldSrcNews[c], localDestBox, localSrcBox, alpha);
        }
    }
};

template<typename GridLayoutT, typename FieldT, typename PhysicalQuantity>
using VecFieldLinearTimeInterpolate
    = TensorFieldLinearTimeInterpolate<1, GridLayoutT, FieldT, PhysicalQuantity>;


} // namespace PHARE::amr

#endif
