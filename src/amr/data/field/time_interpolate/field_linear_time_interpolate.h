#ifndef PHARE_FIELD_LINEAR_TIME_INTERPOLATE_H
#define PHARE_FIELD_LINEAR_TIME_INTERPOLATE_H


// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "data/field/field_data.h"
#include "data/field/field_geometry.h"

#include <SAMRAI/hier/TimeInterpolateOperator.h>


namespace PHARE
{
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




    virtual void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& where,
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


        auto localBox
            = AMRToLocal(static_cast<std::add_const_t<decltype(finalBox)>>(finalBox), ghostBox);


        if constexpr (dim == 1)
        {
            auto iStart = localBox.lower(dirX);
            auto iEnd   = localBox.upper(dirX);

            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                fieldDest(ix) = (1. - alpha) * fieldSrcOld(ix) + alpha * fieldSrcNew(ix);
            }
        }
        else if constexpr (dim == 2)
        {
            auto iStartX = localBox.lower(dirX);
            auto iEndX   = localBox.upper(dirX);

            auto iStartY = localBox.lower(dirY);
            auto iEndY   = localBox.upper(dirY);

            for (auto ix = iStartX; ix <= iEndX; ++ix)
            {
                for (auto iy = iStartY; iy <= iEndY; ++iy)
                {
                    fieldDest(ix, iy)
                        = (1. - alpha) * fieldSrcOld(ix, iy) + alpha * fieldSrcNew(ix, iy);
                }
            }
        }
        else if constexpr (dim == 3)
        {
            auto iStartX = localBox.lower(dirX);
            auto iEndX   = localBox.upper(dirX);

            auto iStartY = localBox.lower(dirY);
            auto iEndY   = localBox.upper(dirY);

            auto iStartZ = localBox.lower(dirZ);
            auto iEndZ   = localBox.upper(dirZ);

            for (auto ix = iStartX; ix <= iEndX; ++ix)
            {
                for (auto iy = iStartY; iy <= iEndY; ++iy)
                {
                    for (auto iz = iStartZ; iz <= iEndZ; ++iz)
                    {
                        fieldDest(ix, iy, iz) = (1. - alpha) * fieldSrcOld(ix, iy, iz)
                                                + alpha * fieldSrcNew(ix, iy, iz);
                    }
                }
            }
        }
        else
        {
            static_assert(dim >= 0 && dim <= 3);
        }
    }




private:
    static std::size_t constexpr dim = GridLayoutImpl::dimension;

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutT, FieldT>;
};



} // namespace PHARE

#endif
