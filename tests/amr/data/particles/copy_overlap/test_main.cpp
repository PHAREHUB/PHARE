#include "data/particles/particles_data.h"
#include "data/particles/particles_data_factory.h"
#include "data/particles/particles_variable.h"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
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

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{dimension, 0},
                                 SAMRAI::hier::Index{dimension, 5}, blockId};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{dimension, 10},
                                   SAMRAI::hier::Index{dimension, 15}, blockId};

    SAMRAI::hier::IntVector ghost{dimension, 1};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    SAMRAI::hier::Patch destPatch{destDomain, patchDescriptor};
    SAMRAI::hier::Patch sourcePatch{sourceDomain, patchDescriptor};

    ParticlesData<1> destPdat{destDomain, ghost};
    ParticlesData<1> sourcePdat{sourceDomain, ghost};


    std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};


    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{sourcePdat.getGhostBox()};
    SAMRAI::hier::Box fillBox{destPdat.getGhostBox()};

    bool overwriteInterior{true};

    // we want the last cell of source (cell AMR 15) to be on top of the
    // ghost cell of dest (cell AMR -1)
    // so the formula is dest.lower() - source.upper() - 1
    SAMRAI::hier::Transformation transformation{
        (destDomain.lower() - sourceDomain.upper()
         - SAMRAI::hier::Index{SAMRAI::hier::IntVector::getOne(dimension)})};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
            *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};
};




TEST_F(AParticlesData1D, PreserveVelocityWhenCopyingWithPeriodics)
{
    Particle<1> particle;

    particle.weight = 1.0;
    particle.charge = 1.0;
    particle.iCell  = {{6}};
    particle.v      = {{1.0, 1.0, 1.0}};

    sourcePdat.interior.push_back(particle);

    destPdat.copy(sourcePdat, *cellOverlap);

    // ghost start of pach1 is 9, so iCell 0 and 5 correspond to 15.
    // cell 15 correspond to ghost of patch0

    ASSERT_THAT(destPdat.ghost.size(), Eq(1));
    ASSERT_THAT(destPdat.ghost[0].v, Eq(particle.v));
}




TEST_F(AParticlesData1D, ShiftTheiCellWhenCopyingWithPeriodics)
{
    ParticlesData<1> pDat0{destDomain, ghost};
    ParticlesData<1> pDat1{sourceDomain, ghost};

    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {{1.0, 1.0, 1.0}};

    sourcePdat.interior.push_back(particle1);

    destPdat.copy(sourcePdat, *cellOverlap);

    // patch0 start at 0 , patch1 start at 10
    // with periodics condition, we have -1 equivalent to 15
    std::array<int, 1> expectediCell{{0}};


    ASSERT_THAT(destPdat.ghost.size(), Eq(1));
    ASSERT_THAT(destPdat.ghost[0].iCell, Eq(expectediCell));
}




TEST_F(AParticlesData1D, CopyBorderSourceParticlesIntoDestGhostWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;
    particle1.iCell  = {{6}};
    particle1.v      = {{1.0, 1.0, 1.0}};


    sourcePdat.ghost.push_back(particle1);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghost.size(), Eq(1));
}




TEST_F(AParticlesData1D, CopyGhostSourceParticlesIntoInteriorDestWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;
    particle1.iCell  = {{7}};
    particle1.v      = {{1.0, 1.0, 1.0}};


    sourcePdat.ghost.push_back(particle1);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.interior.size(), Eq(1));
}




TEST_F(AParticlesData1D, PreserveWeightWhenCopyingWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    sourcePdat.interior.push_back(particle1);

    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghost[0].weight, Eq(particle1.weight));
}




TEST_F(AParticlesData1D, PreserveChargeWhenCopyingWithPeriodics)
{
    Particle<1> particle1;

    particle1.weight = 1.0;
    particle1.charge = 1.0;

    particle1.iCell = {{6}};

    particle1.v = {1.0, 1.0, 1.0};

    sourcePdat.interior.push_back(particle1);

    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghost[0].charge, Eq(particle1.charge));
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
