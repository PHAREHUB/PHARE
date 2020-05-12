#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_H
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_H


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "amr/data/field/field_data.h"
#include "amr/data/field/field_geometry.h"

#include <SAMRAI/hier/TimeInterpolateOperator.h>


namespace PHARE
{
namespace amr
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




        virtual void timeInterpolate(SAMRAI::hier::PatchData& destData,
                                     SAMRAI::hier::Box const& where,
                                     SAMRAI::hier::BoxOverlap const& /*overlap*/,
                                     SAMRAI::hier::PatchData const& srcDataOld,
                                     SAMRAI::hier::PatchData const& srcDataNew) const override
        {
            //

            auto& fieldDataDest = dynamic_cast<FieldDataT&>(destData);

            auto const& fieldDataSrcOld = dynamic_cast<FieldDataT const&>(srcDataOld);
            auto const& fieldDataSrcNew = dynamic_cast<FieldDataT const&>(srcDataNew);

            double interpTime = fieldDataDest.getTime();

            double oldTime = fieldDataSrcOld.getTime();
            double newTime = fieldDataSrcNew.getTime();


            double alpha = (interpTime - oldTime) / (newTime - oldTime);

            auto& layout = fieldDataDest.gridLayout;



            auto& fieldDest = fieldDataDest.field;

            auto& fieldSrcOld = fieldDataSrcOld.field;
            auto& fieldSrcNew = fieldDataSrcNew.field;



            auto qty = fieldDest.physicalQuantity();


            auto whereLayout
                = FieldGeometry<GridLayoutT, PhysicalQuantity>::layoutFromBox(where, layout);

            bool const withGhost{true};
            auto interpolateBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                where, qty, whereLayout, !withGhost);

            auto ghostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                fieldDataDest.getBox(), qty, layout, withGhost);

            auto finalBox = interpolateBox * ghostBox;

            auto srcGhostBox = FieldGeometry<GridLayoutT, PhysicalQuantity>::toFieldBox(
                fieldDataSrcNew.getBox(), qty, fieldDataSrcNew.gridLayout, withGhost);

            auto localDestBox
                = AMRToLocal(static_cast<std::add_const_t<decltype(finalBox)>>(finalBox), ghostBox);

            auto localSrcBox = AMRToLocal(
                static_cast<std::add_const_t<decltype(finalBox)>>(finalBox), srcGhostBox);


            if constexpr (dim == 1)
            {
                auto iDestStart = localDestBox.lower(dirX);
                auto iDestEnd   = localDestBox.upper(dirX);

                auto iSrcStart = localSrcBox.lower(dirX);

                for (auto ix = iDestStart, ixSrc = iSrcStart; ix <= iDestEnd; ++ix, ++ixSrc)
                {
                    fieldDest(ix) = (1. - alpha) * fieldSrcOld(ixSrc) + alpha * fieldSrcNew(ixSrc);
                }
            }
            else if constexpr (dim == 2)
            {
                auto iDestStartX = localDestBox.lower(dirX);
                auto iDestEndX   = localDestBox.upper(dirX);
                auto iDestStartY = localDestBox.lower(dirY);
                auto iDestEndY   = localDestBox.upper(dirY);

                auto iSrcStartX = localSrcBox.lower(dirX);
                auto iSrcStartY = localSrcBox.lower(dirY);

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
                auto iStartX = localDestBox.lower(dirX);
                auto iEndX   = localDestBox.upper(dirX);
                auto iStartY = localDestBox.lower(dirY);
                auto iEndY   = localDestBox.upper(dirY);
                auto iStartZ = localDestBox.lower(dirZ);
                auto iEndZ   = localDestBox.upper(dirZ);

                auto iSrcStartX = localSrcBox.lower(dirX);
                auto iSrcStartY = localSrcBox.lower(dirY);
                auto iSrcStartZ = localSrcBox.lower(dirZ);

                for (auto ix = iStartX, ixSrc = iSrcStartX; ix <= iEndX; ++ix, ++ixSrc)
                {
                    for (auto iy = iStartY, iySrc = iSrcStartY; iy <= iEndY; ++iy, ++ixSrc)
                    {
                        for (auto iz = iStartZ, izSrc = iSrcStartZ; iz <= iEndZ; ++iz, ++izSrc)
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

} // namespace amr

} // namespace PHARE

#endif
