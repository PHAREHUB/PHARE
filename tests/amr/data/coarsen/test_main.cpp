#include <SAMRAI/hier/Box.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>


#include "data/coarsening/field_coarsen.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include <fstream>
#include <numeric>


#include "test_basic_hierarchy.h"
#include "test_linear_coarsen.h"
#include "test_weighter_coarsen.h"



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);



    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);

    SAMRAI::tbox::SAMRAIManager::initialize();

    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();


    // Finalize

    SAMRAI::tbox::SAMRAIManager::shutdown();

    SAMRAI::tbox::SAMRAIManager::finalize();

    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
