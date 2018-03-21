#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <tools/resources_manager.h>

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include <string>
#include <utility>
#include <vector>

#include "basic_hierarchy.h"
#include "cell_placeholder.h"
#include "resource_test_1d.h"

std::string inputBase;



int main(int argc, char **argv)
{
    if (argc > 1)
    {
        inputBase = std::string{argv[1]};
    }
    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    ::testing::InitGoogleTest(&argc, argv);

    int testResult = RUN_ALL_TESTS();

    // Finalize
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
