
#include "core/def/phare_mpi.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "amr/data/particles/particles_data.hpp"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using testing::DoubleEq;
using testing::Eq;
using testing::Pointwise;

using namespace PHARE;
using namespace PHARE::core;
using namespace PHARE::amr;

template<typename Particles>
auto particle_fetcher(Particles const& particles, std::size_t const idx = 0)
{
    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, AoSTS))
    {
        for (auto const& tile : particles())
        {
            if (tile().size())
                return tile()[idx];
        }
        throw std::runtime_error("no available particle found");
    }
    else
        return particles[idx];
}

template<std::size_t _dim, auto _layout_mode, auto _alloc_mode>
struct TestParam
{
    static_assert(all_are<LayoutMode>(_layout_mode));
    static_assert(all_are<AllocatorMode>(_alloc_mode));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
};

template<typename Param>
struct AParticlesData
{
    static constexpr auto dim         = Param::dim;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;

    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    using Particle_t = Particle<dim>;

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

    ParticlesData<ParticleArray_t> destData{destDomain, ghost, "name"};
    ParticlesData<ParticleArray_t> sourceData{sourceDomain, ghost, "name"};
    Particle_t particle;

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


// clang-format off
using ParticlesDatas = testing::Types<
    AParticlesData<TestParam<1, LayoutMode::AoSMapped, AllocatorMode::CPU>>,
    AParticlesData<TestParam<2, LayoutMode::AoS, AllocatorMode::CPU>>,
    AParticlesData<TestParam<2, LayoutMode::AoSTS, AllocatorMode::CPU>>,
    AParticlesData<TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>>
>;
// clang-format on

TYPED_TEST_SUITE(CopyTest, ParticlesDatas);



TYPED_TEST(CopyTest, copiesSourceDomainParticleIntoDomainDestForDomainOverlapCells)
{
    static constexpr auto dim = TypeParam::dim;

    // this set of particles is in the domain of the source patchdata
    // and in domain of the destination patchdata
    for (auto iCell = 3; iCell <= 5; ++iCell)
    {
        this->particle.iCell_ = ConstArray<int, dim>(iCell);

        ASSERT_THAT(this->sourceData.domainParticles.size(), Eq(0));
        this->sourceData.domainParticles.push_back(this->particle);
        ASSERT_THAT(this->sourceData.domainParticles.size(), Eq(1));

        ASSERT_THAT(this->destData.domainParticles.size(), Eq(0));
        ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
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

    auto const& domainParticle = particle_fetcher(this->destData.domainParticles);
    EXPECT_THAT(domainParticle.iCell(), Eq(this->particle.iCell()));
    EXPECT_THAT(domainParticle.delta(), Pointwise(DoubleEq(), this->particle.delta()));

    EXPECT_THAT(domainParticle.v(), Pointwise(DoubleEq(), this->particle.v()));
    EXPECT_THAT(domainParticle.weight(), DoubleEq(this->particle.weight()));
    EXPECT_THAT(domainParticle.charge(), DoubleEq(this->particle.charge()));

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata

    // <<<<<<< HEAD
    //     this->sourceData.domainParticles.clear();
    //     this->destData.domainParticles.clear();

    //     EXPECT_THAT(this->sourceData.domainParticles.size(), Eq(0));
    //     EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    //     EXPECT_THAT(this->sourceData.patchGhostParticles.size(), Eq(0));
    //     EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));

    //     auto const newCell   = ConstArray<int, dim>(6);
    //     this->particle.iCell = ConstArray<int, dim>(6);
    //     EXPECT_THAT(newCell, Eq(this->particle.iCell));

    // =======
    //     this->particle.iCell_              = ConstArray<int, dim>(6);
    // >>>>>>> 0cc9019 (...)
    //     this->sourceData.domainParticles.push_back(this->particle);
    //     EXPECT_THAT(this->sourceData.domainParticles.size(), Eq(1));
    //     EXPECT_THAT(this->sourceData.domainParticles[0].iCell, Eq(newCell));

    // <<<<<<< HEAD
    //     this->destData.copy(this->sourceData);
    //     EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));

    //     EXPECT_THAT(this->destData.domainParticles[0].v, Pointwise(DoubleEq(),
    //     this->particle.v)); EXPECT_THAT(this->destData.domainParticles[0].iCell, Eq(newCell));
    //     EXPECT_THAT(this->destData.domainParticles[0].iCell, Eq(this->particle.iCell));
    //     EXPECT_THAT(this->destData.domainParticles[0].delta,
    //                 Pointwise(DoubleEq(), this->particle.delta));
    //     EXPECT_THAT(this->destData.domainParticles[0].weight, DoubleEq(this->particle.weight));
    //     EXPECT_THAT(this->destData.domainParticles[0].charge, DoubleEq(this->particle.charge));
    // =======
    //     auto const& destPatchGhostParticle =
    //     particle_fetcher(this->destData.patchGhostParticles);
    //     EXPECT_THAT(destPatchGhostParticle.iCell(), Eq(this->particle.iCell()));
    //     EXPECT_THAT(destPatchGhostParticle.delta(), Pointwise(DoubleEq(),
    //     this->particle.delta())); EXPECT_THAT(destPatchGhostParticle.v(), Pointwise(DoubleEq(),
    //     this->particle.v())); EXPECT_THAT(destPatchGhostParticle.weight(),
    //     DoubleEq(this->particle.weight())); EXPECT_THAT(destPatchGhostParticle.charge(),
    //     DoubleEq(this->particle.charge()));
    // >>>>>>> 0cc9019 (...)
}


TYPED_TEST(CopyTest, copiesDataWithOverlapNoTransform)
{
    static constexpr auto dim = TypeParam::dim;
    auto dimension            = SAMRAI::tbox::Dimension{dim};

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


    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(3);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    this->particle.iCell_ = ConstArray<int, dim>(6);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));

    this->sourceData.domainParticles.clear();
    this->sourceData.patchGhostParticles.clear();
    this->destData.patchGhostParticles.clear();
    this->destData.domainParticles.clear();
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(6);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));

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

    // particle is in the domain of the source patchdata
    // and in domain of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(7);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    EXPECT_EQ(5, particle_fetcher(this->destData.domainParticles).iCell()[0]);

    this->sourceData.domainParticles.clear();
    this->destData.domainParticles.clear();

    // particle is in the domain of the source patchdata
    // and in last ghost of the destination patchdata
    // and also in the overlap

    this->particle.iCell_ = ConstArray<int, dim>(8);

    this->sourceData.domainParticles.push_back(this->particle);
    this->destData.copy(this->sourceData, overlap);
    // <<<<<<< HEAD
    //     EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
    //     EXPECT_THAT(this->destData.domainParticles.size(), Eq(1));
    //     EXPECT_EQ(6, this->destData.domainParticles[0].iCell[0]);
    // =======
    //     EXPECT_THAT(this->destData.patchGhostParticles.size(), Eq(1));
    //     EXPECT_THAT(this->destData.domainParticles.size(), Eq(0));
    //     EXPECT_EQ(6, particle_fetcher(this->destData.patchGhostParticles).iCell()[0]);
    // >>>>>>> 0cc9019 (...)

    this->sourceData.domainParticles.clear();
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
