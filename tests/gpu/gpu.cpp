
#define PHARE_WITH_GPU

#include "gpu.hpp"

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    KLOG(NON) << "float64 start";
    PHARE::gpu::do_thing<double>("job_1d");
    return 0;
}
