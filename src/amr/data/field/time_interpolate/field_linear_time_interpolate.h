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
template<typename GridLayoutImpl, typename FieldT>
class FieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
public:
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


        bool const withGhost{true};
        auto interpolateBox = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
            where, qty, layout, !withGhost);

        auto ghostBox = FieldGeometry<GridLayoutImpl, PhysicalQuantity>::toFieldBox(
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
    }




private:
    static std::size_t constexpr dim = GridLayoutImpl::dimension;

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutImpl, FieldT>;
};



} // namespace PHARE

#endif
