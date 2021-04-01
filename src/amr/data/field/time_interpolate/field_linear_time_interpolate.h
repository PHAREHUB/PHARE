#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_H
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_H


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "amr/data/field/field_data.h"
#include "amr/data/field/field_geometry.h"

#include <SAMRAI/hier/TimeInterpolateOperator.h>


namespace PHARE::amr
{
using core::dirX;
using core::dirY;
using core::dirZ;

template<typename GridLayoutT, typename FieldT>
class FieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
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
        //

        auto& fieldDataDest = dynamic_cast<FieldDataT&>(destData);

        auto const& fieldDataSrcOld = dynamic_cast<FieldDataT const&>(srcDataOld);
        auto const& fieldDataSrcNew = dynamic_cast<FieldDataT const&>(srcDataNew);

        double const interpTime = fieldDataDest.getTime();

        double const oldTime = fieldDataSrcOld.getTime();
        double const newTime = fieldDataSrcNew.getTime();


        double const alpha = (interpTime - oldTime) / (newTime - oldTime);

        auto const& layout = fieldDataDest.gridLayout;



        auto& fieldDest = fieldDataDest.field;

        auto const& fieldSrcOld = fieldDataSrcOld.field;
        auto const& fieldSrcNew = fieldDataSrcNew.field;



        auto qty = fieldDest.physicalQuantity();


        auto const whereLayout
            = FieldGeometry<GridLayoutT, PhysicalQuantity>::layoutFromBox(where, layout);

        bool const withGhost{true};
        auto const interpolateBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            where, qty, whereLayout, !withGhost);

        auto const ghostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            fieldDataDest.getBox(), qty, layout, withGhost);

        auto const finalBox = interpolateBox * ghostBox;

        auto srcGhostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
            fieldDataSrcNew.getBox(), qty, fieldDataSrcNew.gridLayout, withGhost);

        auto const localDestBox
            = AMRToLocal(static_cast<std::add_const_t<decltype(finalBox)>>(finalBox), ghostBox);

        auto const localSrcBox
            = AMRToLocal(static_cast<std::add_const_t<decltype(finalBox)>>(finalBox), srcGhostBox);


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
                    fieldDest(ix, iy) = (1. - alpha) * fieldSrcOld(ixSrc, iySrc)
                                        + alpha * fieldSrcNew(ixSrc, iySrc);
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
        else
        {
            static_assert(dim > 0 && dim <= 3);
        }
    }




private:
    static std::size_t constexpr dim = GridLayoutImpl::dimension;

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutT, FieldT>;
};

} // namespace PHARE::amr

#endif
