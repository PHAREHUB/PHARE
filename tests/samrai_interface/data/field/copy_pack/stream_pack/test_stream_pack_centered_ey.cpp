#include "test_stream_pack_centered_ey.h"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr_interface;


TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy);


TYPED_TEST_P(AFieldData1DCenteredOnEy, PackStreamLikeANodeData)
{
    auto& destinationLayout = this->param.destinationFieldData->gridLayout;



    int lower = 6;
    int upper = 9;

    if (destinationLayout.interp_order >= 2)
    {
        upper = 15;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, lower},
                               SAMRAI::hier::Index{this->dim, upper}, this->blockId};

    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(this->dim)};

    Setup1DCenteredOnEy<TypeParam> testWithBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEy, PackStreamWithPeriodicsLikeANodeData)
{
    auto& param = this->param;
    SAMRAI::hier::IntVector shift{param.destinationPatch.getBox().lower()
                                  - param.sourcePatch.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};


    SAMRAI::hier::Box srcMask{this->sourceNodeData->getBox()};
    SAMRAI::hier::Box fillMask{this->destinationNodeData->getGhostBox()};

    Setup1DCenteredOnEy<TypeParam> testWithBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEy, PackStreamARegionWithPeriodicsLikeANodeData)
{
    auto& param = this->param;

    SAMRAI::hier::IntVector shift{param.destinationPatch.getBox().lower()
                                  - param.sourcePatch.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};


    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, 15},
                              SAMRAI::hier::Index{this->dim, 20}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, 0}, SAMRAI::hier::Index{this->dim, 3},
                               this->blockId};

    Setup1DCenteredOnEy<TypeParam> testWithBox{*this, srcMask, fillMask, transformation};
}




REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy, PackStreamLikeANodeData,
                           PackStreamWithPeriodicsLikeANodeData,
                           PackStreamARegionWithPeriodicsLikeANodeData);


INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEy,
                              FieldDataTestList);
