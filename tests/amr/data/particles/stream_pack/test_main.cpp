
#include "phare_mpi.hpp"

#include "core/utilities/types.hpp"

#include "amr/data/particles/particles_data.hpp"

#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "tests/amr/amr.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>

using testing::Eq;

using namespace PHARE;
using namespace PHARE::core;
using namespace PHARE::amr;


template<std::size_t _dim, auto _layout_mode, auto _alloc_mode>
struct TestParam
{
    static_assert(all_are<LayoutMode>(_layout_mode));
    static_assert(all_are<AllocatorMode>(_alloc_mode));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = _layout_mode;
    auto constexpr static alloc_mode  = _alloc_mode;
};


template<typename Particles>
auto particle_fetcher(Particles const& particles, std::size_t const idx = 0)
{
    using enum LayoutMode;
    if constexpr (any_in(Particles::layout_mode, AoSPCTS))
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


template<typename Param>
struct AParticlesData
{
    static constexpr auto dim         = Param::dim;
    static constexpr auto dimension   = Param::dim;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;

    static constexpr auto ghosts = 1;

    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    using Particle_t = Particle<dim>;


    SAMRAI::tbox::Dimension amr_dimension{dim};
    SAMRAI::hier::BlockId blockId{0};

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{amr_dimension, 0},
                                 SAMRAI::hier::Index{amr_dimension, 5}, blockId};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{amr_dimension, 10},
                                   SAMRAI::hier::Index{amr_dimension, 15}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(amr_dimension)};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    SAMRAI::hier::Patch destPatch{destDomain, patchDescriptor};
    SAMRAI::hier::Patch sourcePatch{sourceDomain, patchDescriptor};

    ParticlesData<ParticleArray_t> destData{destDomain, ghost, "name"};
    ParticlesData<ParticleArray_t> sourceData{sourceDomain, ghost, "name"};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{sourceData.getGhostBox()};
    SAMRAI::hier::Box fillBox{destData.getGhostBox()};

    bool overwriteInterior{true};

    SAMRAI::hier::Index oneIndex{SAMRAI::hier::IntVector::getOne(amr_dimension)};
    SAMRAI::hier::Transformation transformation{
        SAMRAI::hier::IntVector{destDomain.lower() - sourceDomain.upper() - oneIndex}};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
            *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};


    Particle_t particle;


    AParticlesData(std::array<int, dim> const& iCell)
    {
        particle.weight_ = 1.0;
        particle.v_      = {1.0, 1.0, 1.0};
        particle.iCell() = iCell;
        sourceData.domainParticles.push_back(particle);
    }
};



template<typename ParticlesData>
struct StreamPackTest : public ::testing::Test
{
};

// clang-format off
using ParticlesDatas = testing::Types<
    AParticlesData<TestParam<1, LayoutMode::AoSMapped, AllocatorMode::CPU>>,
    AParticlesData<TestParam<2, LayoutMode::AoS, AllocatorMode::CPU>>,
    AParticlesData<TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>>
>;
// clang-format on


TYPED_TEST_SUITE(StreamPackTest, ParticlesDatas, );

TYPED_TEST(StreamPackTest, PreserveVelocityWhenPackStreamWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param{ConstArray<int, dim>(15)};
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;


    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(particle_fetcher(destData.domainParticles).v(), Eq(particle.v()));
}



TYPED_TEST(StreamPackTest, ShiftTheiCellWhenPackStreamWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param{ConstArray<int, dim>(15)};
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;


    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    // patch0 start at 0 , patch1 start at 10
    // with periodics condition, we have 0 equivalent to 15
    auto expectediCell = ConstArray<int, dim>(-1);

    ASSERT_THAT(destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(particle_fetcher(destData.domainParticles).iCell(), Eq(expectediCell));
}



TYPED_TEST(StreamPackTest, PackInTheCorrectBufferWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param{ConstArray<int, dim>(15)};
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;


    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    auto expectediCell = ConstArray<int, dim>(-1);

    ASSERT_THAT(destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(particle_fetcher(destData.domainParticles).iCell(), Eq(expectediCell));
}




TYPED_TEST(StreamPackTest,
           PreserveParticleAttributesWhenPackingWithPeriodicsFromGhostSrcToDomainDest)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param{ConstArray<int, dim>(16)};
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;


    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    auto expectediCell = ConstArray<int, dim>(0);

    EXPECT_THAT(destData.domainParticles.size(), Eq(1));

    auto const& domainParticle = particle_fetcher(destData.domainParticles);
    EXPECT_THAT(domainParticle.v(), Eq(particle.v()));
    EXPECT_THAT(domainParticle.iCell(), Eq(expectediCell));
    EXPECT_THAT(domainParticle.delta(), Eq(particle.delta()));
    EXPECT_THAT(domainParticle.weight(), Eq(particle.weight()));
}



TYPED_TEST(StreamPackTest,
           PreserveParticleAttributesWhenPackingWithPeriodicsFromDomainSrcToGhostDest)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param{ConstArray<int, dim>(15)};
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    auto const expectediCell = ConstArray<int, dim>(-1);

    ASSERT_THAT(destData.domainParticles.size(), Eq(1));

    auto const& ghostParticle = particle_fetcher(destData.domainParticles);
    EXPECT_THAT(ghostParticle.v(), Eq(particle.v()));
    EXPECT_THAT(ghostParticle.iCell(), Eq(expectediCell));
    EXPECT_THAT(ghostParticle.delta(), Eq(particle.delta()));
    EXPECT_THAT(ghostParticle.weight(), Eq(particle.weight()));
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::test::amr::SamraiLifeCycle life{argc, argv};
    return RUN_ALL_TESTS();
}
