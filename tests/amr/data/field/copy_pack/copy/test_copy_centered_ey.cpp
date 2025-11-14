#include "field_data_test_param.hpp"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;



TYPED_TEST_SUITE_P(AFieldData1DCenteredOnEy);

TYPED_TEST_P(AFieldData1DCenteredOnEy, CopyLikeACellData)
{
    // Here we want to copy data from one patch (source)
    // to another patch (destination)
    // we will use 2 different function to fill them
    // sourceFill for the source patch and destinationFill
    // for the destination

    auto& destinationField = this->param.destinationFieldData->field;
    auto& sourceField      = this->param.sourceFieldData->field;


    // we get a pointer to the data (including the ghost region)
    double* destinationNodeStart = this->destinationNodeData->getPointer();

    // since our data is supposed to match a NodeData, we can use our gridlayout
    // to get the correct index within the cellData
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


    // After correctly initializing our data (note that for fieldData it was already performed upon
    // parameter construction)
    this->param.destinationFieldData->copy(*this->param.sourceFieldData);

    this->destinationNodeData->copy(*(this->sourceNodeData));

    iStart = this->param.destinationFieldData->gridLayout.ghostStartIndex(destinationField,
                                                                          Direction::X);
    iEnd   = this->param.destinationFieldData->gridLayout.ghostEndIndex(destinationField,
                                                                        Direction::X);


    double const* nodeDataStart = this->destinationNodeData->getPointer();

    // finnaly we want to check if we have the same values
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(destinationField(ix), Eq(nodeDataStart[ix]));
    }
}

REGISTER_TYPED_TEST_SUITE_P(AFieldData1DCenteredOnEy, CopyLikeACellData);


INSTANTIATE_TYPED_TEST_SUITE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEy,
                               FieldDataTestList);
