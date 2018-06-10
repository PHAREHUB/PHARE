#include "data/particles/particles_data.h"
#include "data/particles/particles_data_factory.h"
#include "data/particles/particles_variable.h"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

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
    Particle<1> particle;

    AParticlesData1D()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {{1.0, 1.0, 1.0}};
    }
};




TEST_F(AParticlesData1D, PreserveVelocityWhenCopyingWithPeriodics)
{
    particle.iCell = {{6}};
    sourcePdat.domainParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghostParticles.size(), Eq(1));
    ASSERT_THAT(destPdat.ghostParticles[0].v, Eq(particle.v));
}




TEST_F(AParticlesData1D, ShiftTheiCellWhenCopyingWithPeriodics)
{
    particle.iCell = {{6}};
    sourcePdat.domainParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    std::array<int, 1> expectediCell{{0}};

    ASSERT_THAT(destPdat.ghostParticles.size(), Eq(1));
    ASSERT_THAT(destPdat.ghostParticles[0].iCell, Eq(expectediCell));
}




TEST_F(AParticlesData1D, CopyBorderSourceParticlesIntoDestGhostWithPeriodics)
{
    particle.iCell = {{6}};

    sourcePdat.ghostParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghostParticles.size(), Eq(1));
}




TEST_F(AParticlesData1D, CopyGhostSourceParticlesIntoInteriorDestWithPeriodics)
{
    particle.iCell = {{7}};
    sourcePdat.ghostParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.domainParticles.size(), Eq(1));
}




TEST_F(AParticlesData1D, PreserveWeightWhenCopyingWithPeriodics)
{
    particle.iCell = {{6}};
    sourcePdat.domainParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghostParticles[0].weight, Eq(particle.weight));
}




TEST_F(AParticlesData1D, PreserveChargeWhenCopyingWithPeriodics)
{
    particle.iCell = {{6}};
    sourcePdat.domainParticles.push_back(particle);
    destPdat.copy(sourcePdat, *cellOverlap);

    ASSERT_THAT(destPdat.ghostParticles[0].charge, Eq(particle.charge));
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
