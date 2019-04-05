#include "data/particles/particles_data.h"
#include "data/particles/particles_data_factory.h"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
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

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{dimension, 0},
                                 SAMRAI::hier::Index{dimension, 5}, blockId};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{dimension, 10},
                                   SAMRAI::hier::Index{dimension, 15}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    SAMRAI::hier::Patch destPatch{destDomain, patchDescriptor};
    SAMRAI::hier::Patch sourcePatch{sourceDomain, patchDescriptor};

    ParticlesData<1> destData{destDomain, ghost};
    ParticlesData<1> sourceData{sourceDomain, ghost};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{sourceData.getGhostBox()};
    SAMRAI::hier::Box fillBox{destData.getGhostBox()};

    bool overwriteInterior{true};

    SAMRAI::hier::Index oneIndex{SAMRAI::hier::IntVector::getOne(dimension)};

    SAMRAI::hier::Transformation transformation{destDomain.lower() - sourceDomain.upper()
                                                - oneIndex};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
            *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};


    Particle<1> particle;


    AParticlesData1D()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {1.0, 1.0, 1.0};
    }
};




TEST_F(AParticlesData1D, PreserveVelocityWhenPackStreamWithPeriodics)
{
    particle.iCell = {{15}};
    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles[0].v, Eq(particle.v));
}




TEST_F(AParticlesData1D, ShiftTheiCellWhenPackStreamWithPeriodics)
{
    particle.iCell = {{15}};

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    // patch0 start at 0 , patch1 start at 10
    // with periodics condition, we have 0 equivalent to 15
    std::array<int, 1> expectediCell{-1};


    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles[0].iCell, Eq(expectediCell));
}




TEST_F(AParticlesData1D, PackInTheCorrectBufferWithPeriodics)
{
    particle.iCell = {{16}};

    sourceData.patchGhostParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    std::array<int, 1> expectediCell{0};

    ASSERT_THAT(destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(destData.domainParticles[0].iCell, Eq(expectediCell));
}




TEST_F(AParticlesData1D, PreserveParticleAttributesWhenPackingWithPeriodicsFromGhostSrcToDomainDest)
{
    particle.iCell = {{16}};

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);


    EXPECT_THAT(destData.domainParticles[0].v, Eq(particle.v));
    EXPECT_THAT(destData.domainParticles[0].iCell[0], Eq(0));
    EXPECT_THAT(destData.domainParticles[0].delta, Eq(particle.delta));
    EXPECT_THAT(destData.domainParticles[0].weight, Eq(particle.weight));
    EXPECT_THAT(destData.domainParticles[0].charge, Eq(particle.charge));
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Ex, particle.Ex);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Ey, particle.Ey);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Ez, particle.Ez);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Bx, particle.Bx);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].By, particle.By);
    EXPECT_DOUBLE_EQ(destData.domainParticles[0].Bz, particle.Bz);
}




TEST_F(AParticlesData1D, PreserveParticleAttributesWhenPackingWithPeriodicsFromDomainSrcToGhostDest)
{
    particle.iCell = {{15}};

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    EXPECT_THAT(destData.patchGhostParticles[0].v, Eq(particle.v));
    EXPECT_THAT(destData.patchGhostParticles[0].iCell[0], Eq(-1));
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
