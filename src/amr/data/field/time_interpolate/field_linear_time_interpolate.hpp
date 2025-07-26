#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_HPP


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "amr/data/field/field_data.hpp"
#include "amr/data/field/field_geometry.hpp"

#include "amr/data/tensorfield/tensor_field_data.hpp"
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/TimeInterpolateOperator.h>
#include <tuple>


namespace PHARE::amr
{

using core::dirX;
using core::dirY;
using core::dirZ;

template<typename Dst>
void linear_time_interpolate(Dst& fieldDest, auto& fieldSrcOld, auto& fieldSrcNew, auto&&... args)
{
    auto static constexpr dim = Dst::dimension;

    auto const& [localDestBox, localSrcBox, alpha] = std::forward_as_tuple(args...);

    if constexpr (dim == 1)
    {
        auto const iDestStartX = localDestBox.lower(dirX);
        auto const iDestEndX   = localDestBox.upper(dirX);

        auto const iSrcStartX = localSrcBox.lower(dirX);

        for (auto ix = iDestStartX, ixSrc = iSrcStartX; ix <= iDestEndX; ++ix, ++ixSrc)
        {
            fieldDest(ix) = (1. - alpha) * fieldSrcOld(ixSrc) + alpha * fieldSrcNew(ixSrc);
        }
    }
    else if constexpr (dim == 2)
    {
        auto const iDestStartX = localDestBox.lower(dirX);
        auto const iDestEndX   = localDestBox.upper(dirX);
        auto const iDestStartY = localDestBox.lower(dirY);
        auto const iDestEndY   = localDestBox.upper(dirY);

        auto const iSrcStartX = localSrcBox.lower(dirX);
        auto const iSrcStartY = localSrcBox.lower(dirY);

        for (auto ix = iDestStartX, ixSrc = iSrcStartX; ix <= iDestEndX; ++ix, ++ixSrc)
        {
            for (auto iy = iDestStartY, iySrc = iSrcStartY; iy <= iDestEndY; ++iy, ++iySrc)
            {
                fieldDest(ix, iy)
                    = (1. - alpha) * fieldSrcOld(ixSrc, iySrc) + alpha * fieldSrcNew(ixSrc, iySrc);
            }
        }
    }
    else if constexpr (dim == 3)
    {
        auto const iDestStartX = localDestBox.lower(dirX);
        auto const iDestEndX   = localDestBox.upper(dirX);
        auto const iDestStartY = localDestBox.lower(dirY);
        auto const iDestEndY   = localDestBox.upper(dirY);
        auto const iDestStartZ = localDestBox.lower(dirZ);
        auto const iDestEndZ   = localDestBox.upper(dirZ);

        auto const iSrcStartX = localSrcBox.lower(dirX);
        auto const iSrcStartY = localSrcBox.lower(dirY);
        auto const iSrcStartZ = localSrcBox.lower(dirZ);

        for (auto ix = iDestStartX, ixSrc = iSrcStartX; ix <= iDestEndX; ++ix, ++ixSrc)
        {
            for (auto iy = iDestStartY, iySrc = iSrcStartY; iy <= iDestEndY; ++iy, ++iySrc)
            {
                for (auto iz = iDestStartZ, izSrc = iSrcStartZ; iz <= iDestEndZ; ++iz, ++izSrc)
                {
                    fieldDest(ix, iy, iz) = (1. - alpha) * fieldSrcOld(ixSrc, iySrc, izSrc)
                                            + alpha * fieldSrcNew(ixSrc, iySrc, izSrc);
                }
            }
        }
    }

    //
}


} // namespace PHARE::amr

namespace PHARE::amr
{

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
