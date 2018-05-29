#include "data/particles/particles_data.h"
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using testing::Eq;

using namespace PHARE;

struct AParticlesData1D : public testing::Test
{
    SAMRAI::tbox::Dimension dimension{1};
    SAMRAI::hier::BlockId blockId{0};

    SAMRAI::hier::Box domain0{SAMRAI::hier::Index{dimension, 1}, SAMRAI::hier::Index{dimension, 5},
                              blockId};

    SAMRAI::hier::Box domain1{SAMRAI::hier::Index{dimension, 3}, SAMRAI::hier::Index{dimension, 10},
                              blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};
};



TEST_F(AParticlesData1D, PreserveVelocityWhenCopying)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{0}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    pDat0.copy(pDat1);

    ASSERT_THAT(pDat0.interior[0].v, Eq(particle1.v));
}

TEST_F(AParticlesData1D, ShiftTheiCellWhenCopying)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{0}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    pDat0.copy(pDat1);

    // patch0 physical start at 1 , patch1 physical start at 3
    // so the origin of patch 1 in patch0 coordinate
    // is 2. Since particle1 start at the origin of patch1,
    // it will be shifted to 2
    std::array<int, 1> expectediCell{2};

    ASSERT_THAT(pDat0.interior[0].iCell, Eq(expectediCell));
}
TEST_F(AParticlesData1D, CopyInTheCorrectBuffer)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    // iCell = 4, originPatch1=2
    // so final particle position = 6
    // patch0 interior is from 1 to 5
    // ghost extend it to 6

    particle1.iCell = {{4}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    pDat0.copy(pDat1);

    ASSERT_THAT(pDat0.ghost.size(), Eq(1));
}
TEST_F(AParticlesData1D, PreserveWeightWhenCopying)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{0}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    pDat0.copy(pDat1);

    ASSERT_THAT(pDat0.interior[0].weight, Eq(particle1.weight));
}
TEST_F(AParticlesData1D, PreserveChargeWhenCopying)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{0}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    pDat0.copy(pDat1);

    ASSERT_THAT(pDat0.interior[0].charge, Eq(particle1.charge));
}

TEST_F(AParticlesData1D, DoNothingForParticleOutOfGhostRegionWhenCopying)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{8}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    pDat0.copy(pDat1);

    EXPECT_THAT(pDat0.ghost.size(), Eq(0));
    EXPECT_THAT(pDat0.interior.size(), Eq(0));
    EXPECT_THAT(pDat0.incoming.size(), Eq(0));
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);



    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);

    SAMRAI::tbox::SAMRAIManager::initialize();

    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();


    // Finalize

    SAMRAI::tbox::SAMRAIManager::shutdown();

    SAMRAI::tbox::SAMRAIManager::finalize();

    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
