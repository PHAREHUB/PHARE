#include "amr/data/particles/particles_data.hpp"
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using testing::DoubleEq;
using testing::Eq;
using testing::Pointwise;

using namespace PHARE::core;
using namespace PHARE::amr;




template<std::size_t dim_, bool SOA = false>
struct AParticlesData
{
    static constexpr auto dim = dim_;
    using ParticleArray_t     = ParticleArray<dim, SOA>;

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

    ParticlesData<ParticleArray_t> destData{destDomain, ghost};
    ParticlesData<ParticleArray_t> sourceData{sourceDomain, ghost};
    Particle<dim> particle;


    AParticlesData()
    {
        particle.weight_ = 1.0;
        particle.charge_ = 1.0;
        particle.v_      = {1.0, 1.0, 1.0};
    }
};

template<typename ParticlesData>
struct CopyTest : public ::testing::Test, public ParticlesData
{
};


using ParticlesDatas
    = testing::Types<AParticlesData<1>, AParticlesData<2>, AParticlesData<3>,
                     AParticlesData<1, true>, AParticlesData<2, true>, AParticlesData<3, true>>;


TYPED_TEST_SUITE(CopyTest, ParticlesDatas);



// TYPED_TEST(CopyTest, copiesSourceGhostParticleIntoDomainForGhostSrcOverDomainDest)
//{
//    static constexpr auto dim = TypeParam::dim;
//
//    // particle is in the first ghost of the source patchdata
//    // and in domain of the destination patchdata
//
//    this->particle.iCell_ = ConstArray<int, dim>(2);
//
//    this->sourceData.patchGhostParticles.push_back(this->particle);
//    this->destData.copy(this->sourceData);
//
//    ASSERT_THAT(this->destData.domainParticles.size(), Eq(1));
//    ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
//}


TYPED_TEST(CopyTest, copiesSourceDomainParticleIntoGhostForDomainSrcOverGhostDest)
{
    static constexpr auto dim = TypeParam::dim;

    // particle is in the domain of the source patchdata
    // and in first ghost of the destination patchdata

    this->particle.iCell_ = ConstArray<int, dim>(6);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(this->destData.domainParticles.size(), Eq(0));
}


TYPED_TEST(CopyTest, copiesSourceDomainParticleIntoDomainDestForDomainOverlapCells)
{
    static constexpr auto dim = TypeParam::dim;

    // this set of particles is in the domain of the source patchdata
    // and in domain of the destination patchdata
    for (auto iCell = 3; iCell <= 5; ++iCell)
    {
        this->particle.iCell_ = ConstArray<int, dim>(iCell);

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


TYPED_TEST(CopyTest, PreservesAllParticleAttributesAfterCopy)
{
    static constexpr auto dim = TypeParam::dim;

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata

    this->particle.iCell_ = ConstArray<int, dim>(3);
    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    auto& destDomainParticles = this->destData.domainParticles;
    EXPECT_THAT(destDomainParticles.iCell(0), Eq(this->particle.iCell()));
    EXPECT_THAT(destDomainParticles.delta(0), Pointwise(DoubleEq(), this->particle.delta()));

    EXPECT_THAT(destDomainParticles.v(0), Pointwise(DoubleEq(), this->particle.v()));
    EXPECT_THAT(destDomainParticles.weight(0), DoubleEq(this->particle.weight()));
    EXPECT_THAT(destDomainParticles.charge(0), DoubleEq(this->particle.charge()));

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata

    this->particle.iCell_ = ConstArray<int, dim>(6);
    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    auto& destPatchGhostParticles = this->destData.patchGhostParticles;
    EXPECT_THAT(destPatchGhostParticles.iCell(0), Eq(this->particle.iCell()));
    EXPECT_THAT(destPatchGhostParticles.delta(0), Pointwise(DoubleEq(), this->particle.delta()));

    EXPECT_THAT(destPatchGhostParticles.v(0), Pointwise(DoubleEq(), this->particle.v()));
    EXPECT_THAT(destPatchGhostParticles.weight(0), DoubleEq(this->particle.weight()));
    EXPECT_THAT(destPatchGhostParticles.charge(0), DoubleEq(this->particle.charge()));
}


TYPED_TEST(CopyTest, copiesDataWithOverlapNoTransform)
{
    static constexpr auto dim = TypeParam::dim;
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

    // this->particle.iCell_ = ConstArray<int, dim>(2);
    //
    //    this->sourceData.patchGhostParticles.push_back(this->particle);
    //    EXPECT_THAT(this->sourceData.domainParticles.size(), Eq(0));
    //    EXPECT_THAT(this->sourceData.patchGhostParticles.size(), Eq(1));
    //
    //    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    //    this->destData.copy(this->sourceData, overlap);
    //    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    //
    //    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    //
    //    this->sourceData.domainParticles.clear();
    //    this->sourceData.patchGhostParticles.clear();
    //    this->destData.patchGhostParticles.clear();
    //    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(3);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    this->particle.iCell_ = {{6}};
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(6);

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

    this->particle.iCell_ = ConstArray<int, dim>(4);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
}



TYPED_TEST(CopyTest, copiesDataWithOverlapWithTransform)
{
    static constexpr auto dim = TypeParam::dim;
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

    //  this->particle.iCell_ = ConstArray<int, dim>(4);
    //
    //    this->sourceData.patchGhostParticles.push_back(this->particle);
    //    this->destData.copy(this->sourceData, overlap);
    //    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    //    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    //    EXPECT_EQ(2, this->destData.domainParticles[0].iCell[0]);
    //
    //    this->sourceData.domainParticles.clear();
    //    this->sourceData.patchGhostParticles.clear();
    //    this->destData.patchGhostParticles.clear();
    //    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(7);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_EQ(5, this->destData.domainParticles.iCell(0)[0]);

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(8);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    EXPECT_EQ(6, this->destData.patchGhostParticles.iCell(0)[0]);

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in the domain of the destination patchdata
    // but not in the overlap... should not be copied in dest

    this->particle.iCell_ = ConstArray<int, dim>(6);

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
