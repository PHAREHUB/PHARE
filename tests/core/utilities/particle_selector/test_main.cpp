
#include <string>

#include "data/particles/particle.h"
#include "utilities/box/box.h"
#include "utilities/particle_selector/particle_selector.h"
#include "utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using PHARE::Box;
using PHARE::Particle;
using PHARE::ParticleSelector;
using PHARE::Point;


TEST(AParticleSelector1D, returnTrueForParticlesInBox)
{
    Box<int, 1> box1{Point<int, 1>{0}, Point<int, 1>{2}};
    Box<int, 1> box2{Point<int, 1>{3}, Point<int, 1>{6}};

    std::vector<Box<int, 1>> boxes;
    boxes.push_back(std::move(box1));
    boxes.push_back(std::move(box2));

    Particle<1> part{0.01, 1, {{4}}, {{0.002f}}, {{1.8, 1.83, 2.28}}};

    auto selector = makeSelector(boxes);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
