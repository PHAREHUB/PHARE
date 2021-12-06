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


template<typename dimType>
struct twoParticlesDataNDTouchingPeriodicBorders : public testing::Test
{
    static constexpr auto dim = dimType{}();

    SAMRAI::tbox::Dimension dimension{dim};
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

    ParticlesData<ParticleArray<dim>> destPdat{destDomain, ghost};
    ParticlesData<ParticleArray<dim>> sourcePdat{sourceDomain, ghost};

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

    Particle<dim> particle;

    twoParticlesDataNDTouchingPeriodicBorders()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {{1.0, 1.0, 1.0}};
    }
};



using WithAllDim = testing::Types<DimConst<1>, DimConst<2>, DimConst<3>>;

TYPED_TEST_SUITE(twoParticlesDataNDTouchingPeriodicBorders, WithAllDim);



TYPED_TEST(twoParticlesDataNDTouchingPeriodicBorders,
           haveATransformationThatPutsUpperSourceCellOnTopOfFirstGhostSourceCell)
{
    EXPECT_EQ(-16, this->transformation.getOffset()[0]);
}


TYPED_TEST(twoParticlesDataNDTouchingPeriodicBorders,
           canCopyUpperSourceParticlesInLowerDestGhostCell)
{
    static constexpr auto dim = TypeParam{}();

    auto leftDestGhostCell = -1;
    auto rightSourceCell   = 15;
    if constexpr (dim == 1)
    {
        this->particle.iCell = {{rightSourceCell}};
    }
    else if constexpr (dim == 2)
    {
        this->particle.iCell = {{rightSourceCell, rightSourceCell}};
    }
    else if constexpr (dim == 3)
    {
        this->particle.iCell = {{rightSourceCell, rightSourceCell, rightSourceCell}};
    }
    this->sourcePdat.domainParticles.push_back(this->particle);
    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));

    EXPECT_THAT(this->destPdat.patchGhostParticles.size(), Eq(1));
    EXPECT_EQ(leftDestGhostCell, this->destPdat.patchGhostParticles[0].iCell[0]);
}


TYPED_TEST(twoParticlesDataNDTouchingPeriodicBorders, preserveParticleAttributesInCopies)
{
    static constexpr auto dim = TypeParam{}();


    this->particle.iCell = ConstArray<int, dim>(15);
    this->sourcePdat.domainParticles.push_back(this->particle);
    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));

    EXPECT_THAT(this->destPdat.patchGhostParticles.size(), Eq(1));
    EXPECT_THAT(this->destPdat.patchGhostParticles[0].v, Eq(this->particle.v));
    EXPECT_THAT(this->destPdat.patchGhostParticles[0].iCell[0], Eq(-1));
    if constexpr (dim > 1)
    {
        EXPECT_THAT(this->destPdat.patchGhostParticles[0].iCell[1], Eq(-1));
    }
    if constexpr (dim > 2)
    {
        EXPECT_THAT(this->destPdat.patchGhostParticles[0].iCell[2], Eq(-1));
    }
    EXPECT_THAT(this->destPdat.patchGhostParticles[0].delta, Eq(this->particle.delta));
    EXPECT_THAT(this->destPdat.patchGhostParticles[0].weight, Eq(this->particle.weight));
    EXPECT_THAT(this->destPdat.patchGhostParticles[0].charge, Eq(this->particle.charge));
}


TYPED_TEST(twoParticlesDataNDTouchingPeriodicBorders,
           CopyGhostSourceParticlesIntoInteriorDestWithPeriodics)
{
    static constexpr auto dim = TypeParam{}();

    auto rightSourceGhostCell = 16;
    auto leftDestCell         = 0;

    this->particle.iCell = ConstArray<int, dim>(rightSourceGhostCell);
    this->sourcePdat.patchGhostParticles.push_back(this->particle);
    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));

    EXPECT_THAT(this->destPdat.domainParticles.size(), Eq(1));
    EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[0]);
    if constexpr (dim > 1)
    {
        EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[1]);
    }
    if constexpr (dim > 2)
    {
        EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[2]);
    }
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
