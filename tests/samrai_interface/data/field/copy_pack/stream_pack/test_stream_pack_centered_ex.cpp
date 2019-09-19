#include "test_stream_pack_centered_ex.h"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr_interface;



TYPED_TEST_SUITE_P(AFieldData1DCenteredOnEx);


TYPED_TEST_P(AFieldData1DCenteredOnEx, PackStreamLikeACellData)
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

    Setup1DCenteredOnEx<TypeParam> testWithBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEx, PackStreamWithPeriodicsLikeACellData)
{
    auto& destinationPatch = this->patch1d.destinationPatch;
    auto& sourcePatch      = this->patch1d.sourcePatch;




    SAMRAI::hier::Transformation transformation{destinationPatch.getBox().lower()
                                                - sourcePatch.getBox().upper()};
    SAMRAI::hier::Box srcMask{this->sourceCellData->getBox()};
    SAMRAI::hier::Box fillMask{this->destinationCellData->getGhostBox()};

    Setup1DCenteredOnEx<TypeParam> testWithBox{*this, srcMask, fillMask, transformation};
}




TYPED_TEST_P(AFieldData1DCenteredOnEx, PackStreamARegionWithPeriodicsLikeACellData)
{
    auto& destinationPatch = this->patch1d.destinationPatch;
    auto& sourcePatch      = this->patch1d.sourcePatch;


    SAMRAI::hier::Transformation transformation{destinationPatch.getBox().lower()
                                                - sourcePatch.getBox().upper()};


    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, 15},
                              SAMRAI::hier::Index{this->dim, 20}, this->blockId};
    SAMRAI::hier::Box fillMask{SAMRAI::hier::Index{this->dim, 0}, SAMRAI::hier::Index{this->dim, 3},
                               this->blockId};


    Setup1DCenteredOnEx<TypeParam> testWithBox{*this, srcMask, fillMask, transformation};
}



REGISTER_TYPED_TEST_SUITE_P(AFieldData1DCenteredOnEx, PackStreamLikeACellData,
                            PackStreamWithPeriodicsLikeACellData,
                            PackStreamARegionWithPeriodicsLikeACellData);


INSTANTIATE_TYPED_TEST_SUITE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEx,
                               FieldDataTestList);
