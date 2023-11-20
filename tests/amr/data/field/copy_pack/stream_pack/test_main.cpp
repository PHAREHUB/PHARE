
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/tbox/SAMRAIManager.hpp>
#include <SAMRAI/tbox/SAMRAI_MPI.hpp>

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
