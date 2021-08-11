

#include "tests/simulator/per_test.h"

template<typename Float, std::size_t dim = 1, std::size_t interp = 1, std::size_t nbRefineParts = 2>
void do_thing()
{
    using PHARE_TYPES = PHARE::PHARE_Types<dim, interp, nbRefineParts, Float>;
    SimulatorTestParam<dim, interp, nbRefineParts, Float> sim{"job_" + std::to_string(dim) + "d"};
}

int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam(argc, argv);
    // KLOG(NON) << "float64 start";
    do_thing<double>();
    // KLOG(NON) << "float32 start";
    do_thing<float>();
    return 0;
}
