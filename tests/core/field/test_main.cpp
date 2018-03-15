#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <core/data/field/field.h>



using PHARE::Field;

TEST(FieldTest, sizeAtCreation)
{
    Field field;

    EXPECT_EQ(field.size(), 0);
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
