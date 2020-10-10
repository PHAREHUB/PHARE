#include "amr/data/particles/particles_data.h"
#include "amr/data/particles/particles_data_factory.h"
#include "amr/data/particles/particles_variable.h"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;

struct twoParticlesDatasTouchingPeriodicBorders : public testing::Test
{
    SAMRAI::tbox::Dimension dimension{1};
    SAMRAI::hier::BlockId blockId{0};

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{dimension, 0},
                                 SAMRAI::hier::Index{dimension, 5}, blockId};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{dimension, 10},
                                   SAMRAI::hier::Index{dimension, 15}, blockId};

    SAMRAI::hier::IntVector ghost{dimension, 1};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    SAMRAI::hier::Patch destPatch{destDomain, patchDescriptor};
    SAMRAI::hier::Patch sourcePatch{sourceDomain, patchDescriptor};

    ParticlesData<ParticleArray<1>> destPdat{destDomain, ghost};
    ParticlesData<ParticleArray<1>> sourcePdat{sourceDomain, ghost};


    std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};


    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{sourcePdat.getGhostBox()};
    SAMRAI::hier::Box fillBox{destPdat.getGhostBox()};

    bool overwriteInterior{true};

    SAMRAI::hier::Transformation transformation{destPdat.getGhostBox().lower()
                                                - sourceDomain.upper()};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
            *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};


    Particle<1> particle;

    twoParticlesDatasTouchingPeriodicBorders()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {{1.0, 1.0, 1.0}};
    }
};


TEST_F(twoParticlesDatasTouchingPeriodicBorders,
       haveATransformationThatPutsUpperSourceCellOnTopOfFirstGhostSourceCell)
{
    EXPECT_EQ(-16, transformation.getOffset()[0]);
}



TEST_F(twoParticlesDatasTouchingPeriodicBorders, canCopyUpperSourceParticlesInLowerDestGhostCell)
{
    auto leftDestGhostCell = -1;
    auto upperSourceCell   = 15;
    particle.iCell         = {{upperSourceCell}};
    sourcePdat.domainParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    EXPECT_THAT(destPdat.patchGhostParticles.size(), Eq(1));
    EXPECT_EQ(leftDestGhostCell, destPdat.patchGhostParticles[0].iCell[0]);
}




TEST_F(twoParticlesDatasTouchingPeriodicBorders, preserveParticleAttributesInCopies)
{
    particle.iCell = {{15}};
    sourcePdat.domainParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    EXPECT_THAT(destPdat.patchGhostParticles.size(), Eq(1));

    EXPECT_THAT(destPdat.patchGhostParticles[0].v, Eq(particle.v));
    // EXPECT_THAT(destPdat.ghostParticles[0].iCell, Eq(-1));
    EXPECT_THAT(destPdat.patchGhostParticles[0].delta, Eq(particle.delta));
    EXPECT_THAT(destPdat.patchGhostParticles[0].weight, Eq(particle.weight));
    EXPECT_THAT(destPdat.patchGhostParticles[0].charge, Eq(particle.charge));
    EXPECT_DOUBLE_EQ(destPdat.patchGhostParticles[0].Ex, particle.Ex);
    EXPECT_DOUBLE_EQ(destPdat.patchGhostParticles[0].Ey, particle.Ey);
    EXPECT_DOUBLE_EQ(destPdat.patchGhostParticles[0].Ez, particle.Ez);
    EXPECT_DOUBLE_EQ(destPdat.patchGhostParticles[0].Bx, particle.Bx);
    EXPECT_DOUBLE_EQ(destPdat.patchGhostParticles[0].By, particle.By);
    EXPECT_DOUBLE_EQ(destPdat.patchGhostParticles[0].Bz, particle.Bz);
}



TEST_F(twoParticlesDatasTouchingPeriodicBorders,
       CopyGhostSourceParticlesIntoInteriorDestWithPeriodics)
{
    auto upperSourceGhostCell = 16;
    auto lowerDestCell        = 0;

    particle.iCell = {{upperSourceGhostCell}};
    sourcePdat.patchGhostParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    EXPECT_THAT(destPdat.domainParticles.size(), Eq(1));
    EXPECT_EQ(lowerDestCell, destPdat.domainParticles[0].iCell[0]);
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
