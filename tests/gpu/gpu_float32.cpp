
#define PHARE_WITH_GPU

#include "gpu.hpp"

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    KLOG(NON) << "float32 start";
    PHARE::gpu::do_thing<float>("job_1d_float32");
    return 0;
}
