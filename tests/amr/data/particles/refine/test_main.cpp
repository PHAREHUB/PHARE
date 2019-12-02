#include "amr/data/particles/particles_data.h"
#include "amr/data/particles/particles_variable.h"
#include "amr/data/particles/refine/particles_data_split.h"
#include "amr/data/particles/refine/split.h"
#include "test_basic_hierarchy.h"
#include "amr/resources_manager/amr_utils.h"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>



#include "gmock/gmock.h"
#include "gtest/gtest.h"


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
