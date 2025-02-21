#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/Box.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>


#include "amr/data/field/coarsening/field_coarsen_index_weight.hpp"

#include "gtest/gtest.h"


#include "test_weighter_coarsen.hpp"



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    int testResult = RUN_ALL_TESTS();
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
    return testResult;
}
