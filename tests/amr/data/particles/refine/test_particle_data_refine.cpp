#include "core/def/phare_mpi.hpp"

#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_variable.hpp"
#include "amr/data/particles/refine/particles_data_split.hpp"
#include "amr/data/particles/refine/split.hpp"
#include "amr/resources_manager/amr_utils.hpp"

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
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
    return testResult;
}
