
#include <string>

#include "core/data/particles/particle.h"
#include "core/utilities/box/box.h"
#include "core/utilities/particle_selector/particle_selector.h"
#include "core/utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using PHARE::core::Box;
using PHARE::core::Particle;
using PHARE::core::ParticleSelectorT;
using PHARE::core::Point;


TEST(AParticleSelector1D, returnTrueForParticlesInBox)
{
    Box<int, 1> box1{Point<int, 1>{0}, Point<int, 1>{2}};
    Box<int, 1> box2{Point<int, 1>{3}, Point<int, 1>{6}};

    std::vector<Box<int, 1>> boxes;
    boxes.push_back(std::move(box1));
    boxes.push_back(std::move(box2));

    [[maybe_unused]] Particle<1> part{0.01, 1, {{4}}, {{0.002f}}, {{1.8, 1.83, 2.28}}};

    [[maybe_unused]] auto selector = makeSelector(boxes);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
