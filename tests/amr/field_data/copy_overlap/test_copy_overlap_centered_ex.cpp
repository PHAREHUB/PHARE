#include "test_copy_overlap_centered_ex.h"

namespace PHARE
{
TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx);



TYPED_TEST_P(AFieldData1DCenteredOnEx, CopyWithOverlapLikeACellData)
{
    int lower = 6;
    int upper = 11;

    auto const& destinationLayout = this->param.destinationFieldData->gridLayout;

    // in case we are at in interp order > 1 we want to see ghost region
    if (destinationLayout.interp_order >= 2)
    {
        upper = 12;
    }

    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(this->dim)};


    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, lower},
                               SAMRAI::hier::Index{this->dim, upper}, this->blockId};

    Setup1DCenteredOnEx<TypeParam> testWithThisBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEx, CopyWithPeriodicsLikeACellData)
{
    SAMRAI::hier::IntVector shift{this->param.destinationPatch.getBox().lower()
                                  - this->param.sourcePatch.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};

    SAMRAI::hier::Box srcMask{this->sourceCellData->getBox()};
    SAMRAI::hier::Box fillMask{this->destinationCellData->getGhostBox()};

    Setup1DCenteredOnEx<TypeParam> testWithThisBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEx, CopyOnARegionWithPeriodicsLikeACellData)
{
    SAMRAI::hier::IntVector shift{this->param.destinationPatch.getBox().lower()
                                  - this->param.sourcePatch.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};



    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, 15},
                              SAMRAI::hier::Index{this->dim, 20}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, 0}, SAMRAI::hier::Index{this->dim, 3},
                               this->blockId};


    Setup1DCenteredOnEx<TypeParam> testWithThisBox{*this, srcMask, fillMask, transformation};
}




REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx, CopyWithOverlapLikeACellData,
                           CopyWithPeriodicsLikeACellData, CopyOnARegionWithPeriodicsLikeACellData);


INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEx,
                              FieldDataTestList);
} // namespace PHARE
