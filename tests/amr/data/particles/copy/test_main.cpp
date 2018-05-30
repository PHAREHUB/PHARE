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

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{dimension, 1},
                                 SAMRAI::hier::Index{dimension, 5}, blockId};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{dimension, 3},
                                   SAMRAI::hier::Index{dimension, 10}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

    ParticlesData<1> destData{destDomain, ghost};
    ParticlesData<1> sourceData{sourceDomain, ghost};
    Particle<1> particle;


    AParticlesData1D()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {1.0, 1.0, 1.0};
    }
};



TEST_F(AParticlesData1D, PreserveVelocityWhenCopying)
{
    particle.iCell = {{0}};
    sourceData.interior.push_back(particle);
    destData.copy(sourceData);

    ASSERT_THAT(destData.interior[0].v, Eq(particle.v));
}




TEST_F(AParticlesData1D, ShiftTheiCellWhenCopying)
{
    particle.iCell = {{0}};
    sourceData.interior.push_back(particle);
    destData.copy(sourceData);

    // patch0 physical start at 1 , patch1 physical start at 3
    // so the origin of patch 1 in patch0 coordinate
    // is 2. Since particle1 start at the origin of patch1,
    // it will be shifted to 2
    std::array<int, 1> expectediCell{2};

    ASSERT_THAT(destData.interior[0].iCell, Eq(expectediCell));
}




TEST_F(AParticlesData1D, CopyInTheCorrectBuffer)
{
    // iCell = 4, originPatch1=2
    // so final particle position = 6
    // patch0 interior is from 1 to 5
    // ghost extend it to 6

    particle.iCell = {{4}};
    sourceData.interior.push_back(particle);
    destData.copy(sourceData);

    ASSERT_THAT(destData.ghost.size(), Eq(1));
}




TEST_F(AParticlesData1D, PreserveWeightWhenCopying)
{
    particle.iCell = {{0}};
    sourceData.interior.push_back(particle);
    destData.copy(sourceData);

    ASSERT_THAT(destData.interior[0].weight, Eq(particle.weight));
}




TEST_F(AParticlesData1D, PreserveChargeWhenCopying)
{
    particle.iCell = {{0}};
    sourceData.interior.push_back(particle);
    destData.copy(sourceData);

    ASSERT_THAT(destData.interior[0].charge, Eq(particle.charge));
}




TEST_F(AParticlesData1D, DoesNothingForParticleOutOfGhostRegionWhenCopying)
{
    particle.iCell = {{8}};
    sourceData.interior.push_back(particle);
    destData.copy(sourceData);

    EXPECT_THAT(destData.ghost.size(), Eq(0));
    EXPECT_THAT(destData.interior.size(), Eq(0));
    EXPECT_THAT(destData.coarseToFine.size(), Eq(0));
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
