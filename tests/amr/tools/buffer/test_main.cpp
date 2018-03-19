#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <amr/tools/buffer.h>



using PHARE::Buffer;

TEST(BufferTest, createBuffer)
{
    Buffer buffer;
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
