
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/data/particles/particle.h"
#include "core/data/particles/particle_array.h"
#include "core/data/particles/particle_utilities.h"
#include "core/utilities/box/box.h"
#include "core/utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include <string>


using namespace PHARE::core;

class AParticle : public ::testing::Test
{
protected:
    Particle<3> part;

public:
    AParticle()
        : part{0.01, 1, {{43, 75, 92}}, {{0.002f, 0.2f, 0.8f}}, {{1.8, 1.83, 2.28}}}
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
    EXPECT_EQ(43, part.iCell[0]);
    EXPECT_EQ(75, part.iCell[1]);
    EXPECT_EQ(92, part.iCell[2]);
}


TEST_F(AParticle, CanBeReducedToAnIntegerPoint)
{
    auto p = cellAsPoint(part);
    EXPECT_EQ(43, p[0]);
    EXPECT_EQ(75, p[1]);
    EXPECT_EQ(92, p[2]);
}



TEST_F(AParticle, CanBeReducedToAnAbsolutePositionPoint)
{
    Point<double, 3> origin;
    std::array<double, 3> meshSize{{0.2, 0.05, 0.4}};
    std::array<std::uint32_t, 3> nbrCells{{20, 30, 40}};
    GridLayout<GridLayoutImplYee<3, 1>> layout{meshSize, nbrCells, origin,
                                               Box{Point{40, 60, 80}, Point{59, 89, 119}}};

    auto iCell            = layout.AMRToLocal(part.iCell);
    auto p                = positionAsPoint(part, layout);
    auto startIndexes     = layout.physicalStartIndex(QtyCentering::primal);
    auto expectedPosition = Point<double, 3>{};
    for (auto i = 0u; i < 3; ++i)
    {
        expectedPosition[i]
            = origin[i] + meshSize[i] * (iCell[i] - startIndexes[i] + part.delta[i]),
            p[i];
        EXPECT_DOUBLE_EQ(expectedPosition[i], p[i]);
    }
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
