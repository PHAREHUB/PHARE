#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <core/data/particles/particle.h>
#include <core/data/particles/particle_array.h>
#include <string>

using PHARE::Particle;

class ParticleTest : public ::testing::Test
{
protected:
    Particle<3> part;

public:
    ParticleTest()
        : part{0.01, 1, {{12, 24, 36}}, {{0.002f, 0.2f, 0.8f}}, {{1.8, 1.83, 2.28}}}
    {
    }
};


TEST_F(ParticleTest, ParticleWeightIsWellInitialized)
{
    EXPECT_DOUBLE_EQ(0.01, part.weight);
}

TEST_F(ParticleTest, ParticleChargeIsInitiliazedOK)
{
    EXPECT_DOUBLE_EQ(1., part.charge);
}

TEST_F(ParticleTest, ParticleFieldsAreInitializedToZero)
{
    EXPECT_DOUBLE_EQ(0.0, part.Ex);
    EXPECT_DOUBLE_EQ(0.0, part.Ey);
    EXPECT_DOUBLE_EQ(0.0, part.Ez);

    EXPECT_DOUBLE_EQ(0.0, part.Bx);
    EXPECT_DOUBLE_EQ(0.0, part.By);
    EXPECT_DOUBLE_EQ(0.0, part.Bz);
}


TEST_F(ParticleTest, ParticleVelocityIsInitializedOk)
{
    EXPECT_DOUBLE_EQ(1.8, part.v[0]);
    EXPECT_DOUBLE_EQ(1.83, part.v[1]);
    EXPECT_DOUBLE_EQ(2.28, part.v[2]);
}

TEST_F(ParticleTest, ParticleDeltaIsInitializedOk)
{
    EXPECT_FLOAT_EQ(0.002f, part.delta[0]);
    EXPECT_FLOAT_EQ(0.2f, part.delta[1]);
    EXPECT_FLOAT_EQ(0.8f, part.delta[2]);
}


TEST_F(ParticleTest, ParticleCellIsInitializedOK)
{
    EXPECT_EQ(12, part.iCell[0]);
    EXPECT_EQ(24, part.iCell[1]);
    EXPECT_EQ(36, part.iCell[2]);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
