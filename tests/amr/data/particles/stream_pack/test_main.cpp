#include "data/particles/particles_data.h"
#include "data/particles/particles_data_factory.h"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
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

    SAMRAI::hier::Box domain0{SAMRAI::hier::Index{dimension, 0}, SAMRAI::hier::Index{dimension, 5},
                              blockId};

    SAMRAI::hier::Box domain1{SAMRAI::hier::Index{dimension, 10},
                              SAMRAI::hier::Index{dimension, 15}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    SAMRAI::hier::Patch patch0{domain0, patchDescriptor};
    SAMRAI::hier::Patch patch1{domain1, patchDescriptor};

    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> particles0Geom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(patch0.getBox(), ghost)};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> particles1Geom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(patch1.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{pDat1.getGhostBox()};
    SAMRAI::hier::Box fillBox{pDat0.getGhostBox()};

    bool overwriteInterior{true};

    SAMRAI::hier::Index oneIndex{SAMRAI::hier::IntVector::getOne(dimension)};

    SAMRAI::hier::Transformation transformation{domain0.lower() - domain1.upper() - oneIndex};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(particles0Geom->calculateOverlap(
            *particles1Geom, srcMask, fillBox, overwriteInterior, transformation))};
};


TEST_F(AParticlesData1D, PreserveVelocityWhenPackStreamWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);


    SAMRAI::tbox::MessageStream particlesWriteStream;

    pDat1.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    pDat0.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(pDat0.ghost.size(), Eq(1));
    ASSERT_THAT(pDat0.ghost[0].v, Eq(particle1.v));
}

TEST_F(AParticlesData1D, ShiftTheiCellWhenPackStreamWithPeriodics)
{
    ParticlesData<1> pDat0{domain0, ghost};
    ParticlesData<1> pDat1{domain1, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    pDat1.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    pDat0.unpackStream(particlesReadStream, *cellOverlap);

    // patch0 start at 0 , patch1 start at 10
    // with periodics condition, we have 0 equivalent to 15
    std::array<int, 1> expectediCell{0};


    ASSERT_THAT(pDat0.ghost.size(), Eq(1));
    ASSERT_THAT(pDat0.ghost[0].iCell, Eq(expectediCell));
}

TEST_F(AParticlesData1D, PackInTheCorrectBufferWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.ghost.push_back(particle1);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    pDat1.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    pDat0.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(pDat0.ghost.size(), Eq(1));
}




TEST_F(AParticlesData1D, PreserveWeightWhenPackingWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    pDat1.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    pDat0.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(pDat0.ghost.size(), Eq(1));
    ASSERT_THAT(pDat0.ghost[0].weight, Eq(particle1.weight));
}




TEST_F(AParticlesData1D, PreserveChargeWhenPackingWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    pDat1.interior.push_back(particle1);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    pDat1.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    pDat0.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(pDat0.ghost.size(), Eq(1));
    ASSERT_THAT(pDat0.ghost[0].charge, Eq(particle1.charge));
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
