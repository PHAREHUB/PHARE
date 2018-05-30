#include "field_data_test_param.h"

using testing::Eq;

namespace PHARE
{
TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy);

TYPED_TEST_P(AFieldData1DCenteredOnEy, CopyLikeACellData)
{
    auto& destinationField = this->param.destinationFieldData->field;
    auto& sourceField      = this->param.sourceFieldData->field;


    double* destinationNodeStart = this->destinationNodeData->getPointer();

    auto iStart = this->param.destinationFieldData->gridLayout.ghostStartIndex(destinationField,
                                                                               Direction::X);
    auto iEnd   = this->param.destinationFieldData->gridLayout.ghostEndIndex(destinationField,
                                                                           Direction::X);

    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        destinationNodeStart[ix] = this->param.destinationFill(ix);
    }

    double* sourceNodeStart = this->sourceNodeData->getPointer();

    iStart = this->param.sourceFieldData->gridLayout.ghostStartIndex(sourceField, Direction::X);
    iEnd   = this->param.sourceFieldData->gridLayout.ghostEndIndex(sourceField, Direction::X);


    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        sourceNodeStart[ix] = this->param.sourceFill(ix);
    }


    this->param.destinationFieldData->copy(*this->param.sourceFieldData);

    this->destinationNodeData->copy(*(this->sourceNodeData));

    iStart = this->param.destinationFieldData->gridLayout.ghostStartIndex(destinationField,
                                                                          Direction::X);
    iEnd   = this->param.destinationFieldData->gridLayout.ghostEndIndex(destinationField,
                                                                      Direction::X);


    double const* nodeDataStart = this->destinationNodeData->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(destinationField(ix), Eq(nodeDataStart[ix]));
    }
}

REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy, CopyLikeACellData);


INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEy,
                              FieldDataTestList);



} // namespace PHARE
