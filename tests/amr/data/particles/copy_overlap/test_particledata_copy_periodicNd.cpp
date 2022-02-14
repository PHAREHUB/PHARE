#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_data_factory.hpp"
#include "amr/data/particles/particles_variable.hpp"
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


template<std::size_t dim_, bool SOA = false>
struct AParticlesData
{
    static constexpr auto dim = dim_;
    using ParticleArray_t     = ParticleArray<dim, SOA>;

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

    ParticlesData<ParticleArray_t> destPdat{destDomain, ghost};
    ParticlesData<ParticleArray_t> sourcePdat{sourceDomain, ghost};

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

    AParticlesData()
    {
        particle.weight_ = 1.0;
        particle.charge_ = 1.0;
        particle.v_      = {{1.0, 1.0, 1.0}};
    }
};

template<typename ParticlesData>
struct CopyOverlapTest : public ::testing::Test, public ParticlesData
{
};

using ParticlesDatas
    = testing::Types<AParticlesData<1>, AParticlesData<2>, AParticlesData<3>,
                     AParticlesData<1, true>, AParticlesData<2, true>, AParticlesData<3, true>>;



TYPED_TEST_SUITE(CopyOverlapTest, ParticlesDatas);


TYPED_TEST(CopyOverlapTest, haveATransformationThatPutsUpperSourceCellOnTopOfFirstGhostSourceCell)
{
    EXPECT_EQ(-16, this->transformation.getOffset()[0]);
}


TYPED_TEST(CopyOverlapTest, canCopyUpperSourceParticlesInLowerDestGhostCell)
{
    static constexpr auto dim = TypeParam::dim;

    auto leftDestGhostCell = -1;
    auto rightSourceCell   = 15;

    this->particle.iCell_ = ConstArray<int, dim>(rightSourceCell);

    this->sourcePdat.domainParticles.push_back(this->particle);
    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));

    EXPECT_THAT(this->destPdat.patchGhostParticles.size(), Eq(1));
    EXPECT_EQ(leftDestGhostCell, this->destPdat.patchGhostParticles.iCell(0)[0]);
}


TYPED_TEST(CopyOverlapTest, preserveParticleAttributesInCopies)
{
    static constexpr auto dim = TypeParam::dim;

    this->particle.iCell_ = ConstArray<int, dim>(15);
    this->sourcePdat.domainParticles.push_back(this->particle);
    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));

    EXPECT_THAT(this->destPdat.patchGhostParticles.size(), Eq(1));

    EXPECT_THAT(this->destPdat.patchGhostParticles.weight(0), Eq(this->particle.weight()));
    EXPECT_THAT(this->destPdat.patchGhostParticles.charge(0), Eq(this->particle.charge()));
    EXPECT_THAT(this->destPdat.patchGhostParticles.iCell(0), Eq(ConstArray<int, dim>(-1)));
    EXPECT_THAT(this->destPdat.patchGhostParticles.delta(0), Eq(this->particle.delta()));
    EXPECT_THAT(this->destPdat.patchGhostParticles.v(0), Eq(this->particle.v()));
}


// TYPED_TEST(CopyOverlapTest,
//           CopyGhostSourceParticlesIntoInteriorDestWithPeriodics)
//{
//    static constexpr auto dim = TypeParam{}();
//
//    auto rightSourceGhostCell = 16;
//    auto leftDestCell         = 0;
//
//    this->particle.iCell_ = ConstArray<int, dim>(rightSourceGhostCell);
//    this->sourcePdat.patchGhostParticles.push_back(this->particle);
//    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));
//
//    EXPECT_THAT(this->destPdat.domainParticles.size(), Eq(1));
//    EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[0]);
//    if constexpr (dim > 1)
//    {
//        EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[1]);
//    }
//    if constexpr (dim > 2)
//    {
//        EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[2]);
//    }
//}


// TYPED_TEST(CopyOverlapTest,
//           CopyGhostSourceParticlesIntoInteriorDestWithPeriodics)
//{
//    static constexpr auto dim = TypeParam{}();
//
//    auto rightSourceGhostCell = 16;
//    auto leftDestCell         = 0;
//
//    this->particle.iCell_ = ConstArray<int, dim>(rightSourceGhostCell);
//    this->sourcePdat.patchGhostParticles.push_back(this->particle);
//    this->destPdat.copy(this->sourcePdat, *(this->cellOverlap));
//
//    EXPECT_THAT(this->destPdat.domainParticles.size(), Eq(1));
//    EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[0]);
//    if constexpr (dim > 1)
//    {
//        EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[1]);
//    }
//    if constexpr (dim > 2)
//    {
//        EXPECT_EQ(leftDestCell, this->destPdat.domainParticles[0].iCell[2]);
//    }
//}



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
