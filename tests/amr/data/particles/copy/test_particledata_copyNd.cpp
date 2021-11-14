#include "amr/data/particles/particles_data.h"
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using testing::DoubleEq;
using testing::Eq;
using testing::Pointwise;

using namespace PHARE::core;
using namespace PHARE::amr;


template<typename dimType>
struct AParticlesDataND : public testing::Test
{
    static constexpr auto dim = dimType{}();

    SAMRAI::tbox::Dimension dimension{dim};
    SAMRAI::hier::BlockId blockId{0};

    // the ctor of SAMRAI boxes make them having the same size in each of their dim
    // ghostsize = 1     domain = "---"     ghost = "==="
    //
    //         0   1   2   3   4   5   6   7   8   9  10
    // src :         |===|---------------------------|===|
    // dst : |===|-------------------|===|

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{dimension, 3},
                                   SAMRAI::hier::Index{dimension, 9}, blockId};

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{dimension, 1},
                                 SAMRAI::hier::Index{dimension, 5}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

    ParticlesData<ParticleArray<dim>> destData{destDomain, ghost};
    ParticlesData<ParticleArray<dim>> sourceData{sourceDomain, ghost};
    Particle<dim> particle;


    AParticlesDataND()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {1.0, 1.0, 1.0};
    }
};



using WithAllDim = testing::Types<DimConst<1>, DimConst<2>, DimConst<3>>;

TYPED_TEST_SUITE(AParticlesDataND, WithAllDim);



TYPED_TEST(AParticlesDataND, copiesSourceGhostParticleIntoDomainForGhostSrcOverDomainDest)
{
    static constexpr auto dim = TypeParam{}();

    // particle is in the first ghost of the source patchdata
    // and in domain of the destination patchdata

    this->particle.iCell = ConstArray<int, dim>(2);

    this->sourceData.patchGhostParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    ASSERT_THAT(this->destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
}


TYPED_TEST(AParticlesDataND, copiesSourceDomainParticleIntoGhostForDomainSrcOverGhostDest)
{
    static constexpr auto dim = TypeParam{}();

    // particle is in the domain of the source patchdata
    // and in first ghost of the destination patchdata

    this->particle.iCell = ConstArray<int, dim>(6);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(this->destData.domainParticles.size(), Eq(0));
}


TYPED_TEST(AParticlesDataND, copiesSourceDomainParticleIntoDomainDestForDomainOverlapCells)
{
    static constexpr auto dim = TypeParam{}();

    // this set of particles is in the domain of the source patchdata
    // and in domain of the destination patchdata
    for (auto iCell = 3; iCell <= 5; ++iCell)
    {
        this->particle.iCell = ConstArray<int, dim>(iCell);

        this->sourceData.domainParticles.push_back(this->particle);
        this->destData.copy(this->sourceData);

        ASSERT_THAT(this->destData.domainParticles.size(), Eq(1));
        ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(0));

        this->sourceData.domainParticles.clear();
        this->sourceData.patchGhostParticles.clear();
        this->destData.patchGhostParticles.clear();
        this->destData.domainParticles.clear();
    }
}


TYPED_TEST(AParticlesDataND, PreservesAllParticleAttributesAfterCopy)
{
    static constexpr auto dim = TypeParam{}();

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata

    this->particle.iCell = ConstArray<int, dim>(3);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    EXPECT_THAT(this->destData.domainParticles[0].v, Pointwise(DoubleEq(), this->particle.v));
    EXPECT_THAT(this->destData.domainParticles[0].iCell, Eq(this->particle.iCell));
    EXPECT_THAT(this->destData.domainParticles[0].delta,
                Pointwise(DoubleEq(), this->particle.delta));
    EXPECT_THAT(this->destData.domainParticles[0].weight, DoubleEq(this->particle.weight));
    EXPECT_THAT(this->destData.domainParticles[0].charge, DoubleEq(this->particle.charge));
    EXPECT_DOUBLE_EQ(this->destData.domainParticles[0].Ex, this->particle.Ex);
    EXPECT_DOUBLE_EQ(this->destData.domainParticles[0].Ey, this->particle.Ey);
    EXPECT_DOUBLE_EQ(this->destData.domainParticles[0].Ez, this->particle.Ez);
    EXPECT_DOUBLE_EQ(this->destData.domainParticles[0].Bx, this->particle.Bx);
    EXPECT_DOUBLE_EQ(this->destData.domainParticles[0].By, this->particle.By);
    EXPECT_DOUBLE_EQ(this->destData.domainParticles[0].Bz, this->particle.Bz);

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata

    this->particle.iCell = ConstArray<int, dim>(6);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    EXPECT_THAT(this->destData.patchGhostParticles[0].v, Pointwise(DoubleEq(), this->particle.v));
    EXPECT_THAT(this->destData.patchGhostParticles[0].iCell, Eq(this->particle.iCell));
    EXPECT_THAT(this->destData.patchGhostParticles[0].delta,
                Pointwise(DoubleEq(), this->particle.delta));
    EXPECT_THAT(this->destData.patchGhostParticles[0].weight, DoubleEq(this->particle.weight));
    EXPECT_THAT(this->destData.patchGhostParticles[0].charge, DoubleEq(this->particle.charge));
    EXPECT_DOUBLE_EQ(this->destData.patchGhostParticles[0].Ex, this->particle.Ex);
    EXPECT_DOUBLE_EQ(this->destData.patchGhostParticles[0].Ey, this->particle.Ey);
    EXPECT_DOUBLE_EQ(this->destData.patchGhostParticles[0].Ez, this->particle.Ez);
    EXPECT_DOUBLE_EQ(this->destData.patchGhostParticles[0].Bx, this->particle.Bx);
    EXPECT_DOUBLE_EQ(this->destData.patchGhostParticles[0].By, this->particle.By);
    EXPECT_DOUBLE_EQ(this->destData.patchGhostParticles[0].Bz, this->particle.Bz);
}


TYPED_TEST(AParticlesDataND, copiesDataWithOverlapNoTransform)
{
    static constexpr auto dim = TypeParam{}();
    auto dimension            = SAMRAI::tbox::Dimension{this->dim};

    // now, with an overlap as union of 2 boxes
    //
    //         0   1   2   3   4   5   6   7   8   9  10
    // src :         |===|---------------------------|===|
    // dst : |===|-------------------|===|
    // ovl :         |-------|   |-------|

    SAMRAI::hier::Box box1{SAMRAI::hier::Index{dimension, 2}, SAMRAI::hier::Index{dimension, 3},
                           this->blockId};
    SAMRAI::hier::Box box2{SAMRAI::hier::Index{dimension, 5}, SAMRAI::hier::Index{dimension, 6},
                           this->blockId};

    SAMRAI::hier::BoxContainer container(box1);
    container.push_back(box2);

    SAMRAI::hier::Transformation transfo{SAMRAI::hier::IntVector::getZero(dimension)};
    SAMRAI::pdat::CellOverlap overlap(container, transfo);

    // particle is in the first ghost of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell = ConstArray<int, dim>(2);

    this->sourceData.patchGhostParticles.push_back(this->particle);
    EXPECT_THAT(this->sourceData.domainParticles.size(), Eq(0));
    EXPECT_THAT(this->sourceData.patchGhostParticles.size(), Eq(1));

    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));

    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell = ConstArray<int, dim>(3);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    this->particle.iCell = {{6}};
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata
    // and also in the overlap

    this->particle.iCell = ConstArray<int, dim>(6);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in the domain of the destination patchdata
    // but not in the overlap... should not be copied in dest

    this->particle.iCell = ConstArray<int, dim>(4);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
}



TYPED_TEST(AParticlesDataND, copiesDataWithOverlapWithTransform)
{
    static constexpr auto dim = TypeParam{}();
    auto dimension            = SAMRAI::tbox::Dimension{dim};

    // using the same overlap as in previous test
    SAMRAI::hier::Box box1{SAMRAI::hier::Index{dimension, 2}, SAMRAI::hier::Index{dimension, 3},
                           this->blockId};
    SAMRAI::hier::Box box2{SAMRAI::hier::Index{dimension, 5}, SAMRAI::hier::Index{dimension, 6},
                           this->blockId};

    SAMRAI::hier::BoxContainer container(box1);
    container.push_back(box2);

    // but then we consider an offset of -2 for the overlap (& no rotation)
    SAMRAI::hier::Transformation transfo{SAMRAI::hier::IntVector(dimension, -2)};
    SAMRAI::pdat::CellOverlap overlap(container, transfo);

    // particle is in the first ghost of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell = ConstArray<int, dim>(4);

    this->sourceData.patchGhostParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_EQ(2, this->destData.domainParticles[0].iCell[0]);

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell = ConstArray<int, dim>(7);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_EQ(5, this->destData.domainParticles[0].iCell[0]);

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata
    // and also in the overlap

    this->particle.iCell = ConstArray<int, dim>(8);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    EXPECT_EQ(6, this->destData.patchGhostParticles[0].iCell[0]);

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in the domain of the destination patchdata
    // but not in the overlap... should not be copied in dest

    this->particle.iCell = ConstArray<int, dim>(6);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
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
