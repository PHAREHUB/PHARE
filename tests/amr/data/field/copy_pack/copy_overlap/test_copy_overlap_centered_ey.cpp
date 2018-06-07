#include "test_copy_overlap_centered_ey.h"


using namespace PHARE;



TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy);

TYPED_TEST_P(AFieldData1DCenteredOnEy, CopyWithOverlapLikeANodeData)
{
    int lower = 6;
    int upper = 9;

    auto const& destinationLayout = this->param.destinationFieldData->gridLayout;

    if (destinationLayout.interp_order >= 2)
    {
        upper = 15;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, lower},
                               SAMRAI::hier::Index{this->dim, upper}, this->blockId};


    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(this->dim)};

    Setup1DCenteredOnEy<TypeParam> testWithThisBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEy, CopyWithPeriodicsLikeANodeData)
{
    SAMRAI::hier::IntVector shift{this->param.destinationPatch.getBox().lower()
                                  - this->param.sourcePatch.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};


    SAMRAI::hier::Box srcMask{this->sourceNodeData->getBox()};
    SAMRAI::hier::Box fillMask{this->destinationNodeData->getGhostBox()};

    Setup1DCenteredOnEy<TypeParam> testWithThisBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEy, CopyOnARegionWithPeriodicsLikeANodeData)
{
    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, 15},
                              SAMRAI::hier::Index{this->dim, 20}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, 0}, SAMRAI::hier::Index{this->dim, 3},
                               this->blockId};


    SAMRAI::hier::IntVector shift{this->param.destinationPatch.getBox().lower()
                                  - this->param.sourcePatch.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};



    Setup1DCenteredOnEy<TypeParam> testWithThisBox{*this, srcMask, fillMask, transformation};
}



REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy, CopyWithOverlapLikeANodeData,
                           CopyWithPeriodicsLikeANodeData, CopyOnARegionWithPeriodicsLikeANodeData);


INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEy,
                              FieldDataTestList);
