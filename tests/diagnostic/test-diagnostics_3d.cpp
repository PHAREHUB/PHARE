
#include <mpi.h>

#include "test_diagnostics.ipp"

static std::string const job_file = "job_3d";
static std::string const out_dir  = "phare_outputs/diags_3d/";

TYPED_TEST(Simulator3dTest, fluid)
{
    fluid_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(Simulator3dTest, particles)
{
    particles_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(Simulator3dTest, electromag)
{
    electromag_test(TypeParam{job_file}, out_dir);
}

TYPED_TEST(Simulator3dTest, allFromPython)
{
    allFromPython_test(TypeParam{job_file}, out_dir);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
