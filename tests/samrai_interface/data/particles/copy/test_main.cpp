#include "data/particles/particles_data.h"
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>




#include "gmock/gmock.h"
#include "gtest/gtest.h"


using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr_interface;

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



TEST_F(AParticlesData1D, copiesSourceGhostParticleIntoDomainForGhostSrcOverDomainDest)
{
    // iCell == 2
    // so that the particle is in the first ghost of the source patchdata
    // and in domain of the destination patchdata

    particle.iCell = {{2}};
    sourceData.patchGhostParticles.push_back(particle);
    destData.copy(sourceData);

    ASSERT_THAT(destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(0));
}



TEST_F(AParticlesData1D, copiesSourceDomainParticleIntoGhostForDomainSrcOverGhostDest)
{
    // iCell == 6
    // so that the particle is in the first ghost of the source patchdata
    // and in domain of the destination patchdata

    particle.iCell = {{6}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData);

    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.domainParticles.size(), Eq(0));
}



TEST_F(AParticlesData1D, copiesSourceDomainParticleIntoDomainDestForDomainOverlapCells)
{
    for (auto iCell = 3; iCell <= 5; ++iCell)
    {
        particle.iCell = {{iCell}};
        sourceData.domainParticles.push_back(particle);
        destData.copy(sourceData);

        ASSERT_THAT(destData.domainParticles.size(), Eq(1));
        ASSERT_THAT(destData.patchGhostParticles.size(), Eq(0));

        sourceData.domainParticles.clear();
        sourceData.patchGhostParticles.clear();
        destData.patchGhostParticles.clear();
        destData.domainParticles.clear();
    }
}



TEST_F(AParticlesData1D, PreservesAllParticleAttributesAfterCopy)
{
    particle.iCell = {{3}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData);

    EXPECT_THAT(destData.domainParticles[0].v, Eq(particle.v));
    EXPECT_THAT(destData.domainParticles[0].iCell, Eq(particle.iCell));
    EXPECT_THAT(destData.domainParticles[0].delta, Eq(particle.delta));
    EXPECT_THAT(destData.domainParticles[0].weight, Eq(particle.weight));
    EXPECT_THAT(destData.domainParticles[0].charge, Eq(particle.charge));
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Ex, particle.Ex);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Ey, particle.Ey);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Ez, particle.Ez);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Bx, particle.Bx);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].By, particle.By);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Bz, particle.Bz);


    particle.iCell = {{6}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData);

    EXPECT_THAT(destData.patchGhostParticles[0].v, Eq(particle.v));
    EXPECT_THAT(destData.patchGhostParticles[0].iCell, Eq(particle.iCell));
    EXPECT_THAT(destData.patchGhostParticles[0].delta, Eq(particle.delta));
    EXPECT_THAT(destData.patchGhostParticles[0].weight, Eq(particle.weight));
    EXPECT_THAT(destData.patchGhostParticles[0].charge, Eq(particle.charge));
    EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].Ex, particle.Ex);
    EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].Ey, particle.Ey);
    EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].Ez, particle.Ez);
    EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].Bx, particle.Bx);
    EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].By, particle.By);
    EXPECT_DOUBLE_EQ(destData.patchGhostParticles[0].Bz, particle.Bz);
}




TEST_F(AParticlesData1D, copiesDataWithOverlapNoTransform)
{
    SAMRAI::hier::Box box1{SAMRAI::hier::Index{dimension, 2}, SAMRAI::hier::Index{dimension, 3},
                           blockId};

    SAMRAI::hier::Box box2{SAMRAI::hier::Index{dimension, 4}, SAMRAI::hier::Index{dimension, 6},
                           blockId};

    SAMRAI::hier::BoxContainer container(box1);
    container.push_back(box2);

    SAMRAI::hier::Transformation transfo{
        SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension{1})};

    SAMRAI::pdat::CellOverlap overlap(container, transfo);

    particle.iCell = {{2}};
    sourceData.patchGhostParticles.push_back(particle);
    destData.copy(sourceData, overlap);
    EXPECT_THAT(destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(destData.patchGhostParticles.size(), Eq(0));

    sourceData.domainParticles.clear();
    sourceData.patchGhostParticles.clear();
    destData.patchGhostParticles.clear();
    destData.domainParticles.clear();

    particle.iCell = {{3}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData, overlap);
    EXPECT_THAT(destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(destData.patchGhostParticles.size(), Eq(0));

    sourceData.domainParticles.clear();
    sourceData.patchGhostParticles.clear();
    destData.patchGhostParticles.clear();
    destData.domainParticles.clear();

    particle.iCell = {{6}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData, overlap);
    EXPECT_THAT(destData.patchGhostParticles.size(), Eq(1));
    EXPECT_THAT(destData.domainParticles.size(), Eq(0));
}




TEST_F(AParticlesData1D, copiesDataWithOverlapWithTransform)
{
    SAMRAI::hier::Box box1{SAMRAI::hier::Index{dimension, 2}, SAMRAI::hier::Index{dimension, 3},
                           blockId};

    SAMRAI::hier::Box box2{SAMRAI::hier::Index{dimension, 4}, SAMRAI::hier::Index{dimension, 6},
                           blockId};

    SAMRAI::hier::BoxContainer container(box1);
    container.push_back(box2);

    SAMRAI::hier::Transformation transfo{SAMRAI::hier::IntVector(SAMRAI::tbox::Dimension{1}, -2)};

    SAMRAI::pdat::CellOverlap overlap(container, transfo);

    particle.iCell = {{4}};
    sourceData.patchGhostParticles.push_back(particle);
    destData.copy(sourceData, overlap);
    EXPECT_THAT(destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(destData.patchGhostParticles.size(), Eq(0));
    EXPECT_EQ(2, destData.domainParticles[0].iCell[0]);

    sourceData.domainParticles.clear();
    sourceData.patchGhostParticles.clear();
    destData.patchGhostParticles.clear();
    destData.domainParticles.clear();

    particle.iCell = {{6}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData, overlap);
    EXPECT_THAT(destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(destData.patchGhostParticles.size(), Eq(0));
    EXPECT_EQ(4, destData.domainParticles[0].iCell[0]);

    sourceData.domainParticles.clear();
    sourceData.patchGhostParticles.clear();
    destData.patchGhostParticles.clear();
    destData.domainParticles.clear();

    particle.iCell = {{8}};
    sourceData.domainParticles.push_back(particle);
    destData.copy(sourceData, overlap);
    EXPECT_THAT(destData.patchGhostParticles.size(), Eq(1));
    EXPECT_THAT(destData.domainParticles.size(), Eq(0));
}




int main(int argc, char** argv)
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
