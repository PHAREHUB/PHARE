// comparing std::map and thrust interop

#include <map>

#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>

#include "gtest/gtest.h"

TYPED_TEST(std_map_gpu_tests, thrust_interops)
{
    std::map<int, int> cpu;
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
