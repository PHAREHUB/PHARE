#include "push_bench.hpp"

int main(int /*argc*/, char** /*argv*/)
{
    PHARE::core::bench::write_raw_unsorted_particles_to_file<3>();
}
