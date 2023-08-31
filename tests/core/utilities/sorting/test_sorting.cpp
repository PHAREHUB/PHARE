#include "core/utilities/sorting.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/sorting.hpp"

#include "gtest/gtest.h"


using namespace PHARE::core;


TEST(Sorting, CountingSort)
{
    int const nbrParticles = 1000;
    Box<int, 1> box{{0}, {99}};
    ParticleArray<1> particles(box, nbrParticles);
    CountingSort<ParticleArray<1>, 1> sorter;

    // when allocation is ok find a way to allocate room enough for
    // multiple sorting
    sorter.set_size(nbrParticles);

    // when ready to sort
    auto cellFattener = [](auto const& particle) { return particle.iCell[0]; };
    sorter.sort(particles, cellFattener);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
