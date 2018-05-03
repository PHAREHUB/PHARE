#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <string>

#include "data/particles/particle.h"
#include "data/particles/particle_array.h"
#include "utilities/point/point.h"

using PHARE::cellAsPoint;
using PHARE::Particle;
using PHARE::Point;

class AParticle : public ::testing::Test
{
protected:
    Particle<3> part;

public:
    AParticle()
        : part{0.01, 1, {{12, 24, 36}}, {{0.002f, 0.2f, 0.8f}}, {{1.8, 1.83, 2.28}}}
    {
    }
};


TEST_F(AParticle, ParticleWeightIsWellInitialized)
{
    EXPECT_DOUBLE_EQ(0.01, part.weight);
}

TEST_F(AParticle, ParticleChargeIsInitiliazedOK)
{
    EXPECT_DOUBLE_EQ(1., part.charge);
}

TEST_F(AParticle, ParticleFieldsAreInitializedToZero)
{
    EXPECT_DOUBLE_EQ(0.0, part.Ex);
    EXPECT_DOUBLE_EQ(0.0, part.Ey);
    EXPECT_DOUBLE_EQ(0.0, part.Ez);

    EXPECT_DOUBLE_EQ(0.0, part.Bx);
    EXPECT_DOUBLE_EQ(0.0, part.By);
    EXPECT_DOUBLE_EQ(0.0, part.Bz);
}


TEST_F(AParticle, ParticleVelocityIsInitializedOk)
{
    EXPECT_DOUBLE_EQ(1.8, part.v[0]);
    EXPECT_DOUBLE_EQ(1.83, part.v[1]);
    EXPECT_DOUBLE_EQ(2.28, part.v[2]);
}

TEST_F(AParticle, ParticleDeltaIsInitializedOk)
{
    EXPECT_FLOAT_EQ(0.002f, part.delta[0]);
    EXPECT_FLOAT_EQ(0.2f, part.delta[1]);
    EXPECT_FLOAT_EQ(0.8f, part.delta[2]);
}


TEST_F(AParticle, ParticleCellIsInitializedOK)
{
    EXPECT_EQ(12, part.iCell[0]);
    EXPECT_EQ(24, part.iCell[1]);
    EXPECT_EQ(36, part.iCell[2]);
}


TEST_F(AParticle, CanBeReducedToAnIntegerPoint)
{
    auto p = cellAsPoint(part);
    EXPECT_EQ(12, p[0]);
    EXPECT_EQ(24, p[1]);
    EXPECT_EQ(36, p[2]);
}



TEST_F(AParticle, CanBeReducedToAnAbsolutePositionPoint)
{
    Point<double, 3> origin;
    std::array<double, 3> meshSize{{0.2, 0.05, 0.4}};

    auto p = positionAsPoint(part, meshSize, origin);
    EXPECT_DOUBLE_EQ(origin[0] + meshSize[0] * (part.iCell[0] + part.delta[0]), p[0]);
    EXPECT_DOUBLE_EQ(origin[1] + meshSize[1] * (part.iCell[1] + part.delta[1]), p[1]);
    EXPECT_DOUBLE_EQ(origin[2] + meshSize[2] * (part.iCell[2] + part.delta[2]), p[2]);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
