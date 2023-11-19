
#include "core/def/phare_mpi.hpp"


#include <memory>
#include <cstdint>

#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_data_factory.hpp"
#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchDescriptor.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "core/utilities/types.hpp"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;


template<std::size_t dim>
struct AParticlesData
{
    static constexpr auto dimension = dim;

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

    ParticlesData<ParticleArray<dim>> destData{destDomain, ghost};
    ParticlesData<ParticleArray<dim>> sourceData{sourceDomain, ghost};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(destPatch.getBox(), ghost)};

    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceGeom{
        std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghost)};


    SAMRAI::hier::Box srcMask{sourceData.getGhostBox()};
    SAMRAI::hier::Box fillBox{destData.getGhostBox()};

    bool overwriteInterior{true};

    SAMRAI::hier::Index oneIndex{SAMRAI::hier::IntVector::getOne(amr_dimension)};

    SAMRAI::hier::Transformation transformation{destDomain.lower() - sourceDomain.upper()
                                                - oneIndex};


    std::shared_ptr<SAMRAI::pdat::CellOverlap> cellOverlap{
        std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(destGeom->calculateOverlap(
            *sourceGeom, srcMask, fillBox, overwriteInterior, transformation))};


    typename ParticleArray<dim>::Particle_t particle;


    AParticlesData()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {1.0, 1.0, 1.0};
    }
};


template<typename ParticlesData>
struct StreamPackTest : public ::testing::Test
{
};

using ParticlesDatas = testing::Types<AParticlesData<1>, AParticlesData<2>, AParticlesData<3>>;
TYPED_TEST_SUITE(StreamPackTest, ParticlesDatas);

TYPED_TEST(StreamPackTest, PreserveVelocityWhenPackStreamWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param;
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    particle.iCell = ConstArray<int, dim>(15);
    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles[0].v, Eq(particle.v));
}




TYPED_TEST(StreamPackTest, ShiftTheiCellWhenPackStreamWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param;
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    particle.iCell = ConstArray<int, dim>(15);

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    // patch0 start at 0 , patch1 start at 10
    // with periodics condition, we have 0 equivalent to 15
    auto expectediCell = ConstArray<int, dim>(-1);


    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles[0].iCell, Eq(expectediCell));
}



TYPED_TEST(StreamPackTest, PackInTheCorrectBufferWithPeriodics)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param;
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    particle.iCell = ConstArray<int, dim>(15);

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    auto expectediCell = ConstArray<int, dim>(-1);

    ASSERT_THAT(destData.patchGhostParticles.size(), Eq(1));
    ASSERT_THAT(destData.patchGhostParticles[0].iCell, Eq(expectediCell));
}



TYPED_TEST(StreamPackTest,
           PreserveParticleAttributesWhenPackingWithPeriodicsFromGhostSrcToDomainDest)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param;
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    particle.iCell = ConstArray<int, dim>(16);

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    auto expectediCell = ConstArray<int, dim>(0);

    EXPECT_THAT(destData.domainParticles[0].v, Eq(particle.v));
    EXPECT_THAT(destData.domainParticles[0].iCell, Eq(expectediCell));
    EXPECT_THAT(destData.domainParticles[0].delta, Eq(particle.delta));
    EXPECT_THAT(destData.domainParticles[0].weight, Eq(particle.weight));
    EXPECT_THAT(destData.domainParticles[0].charge, Eq(particle.charge));
}



TYPED_TEST(StreamPackTest,
           PreserveParticleAttributesWhenPackingWithPeriodicsFromDomainSrcToGhostDest)
{
    using ParticlesData = TypeParam;
    constexpr auto dim  = ParticlesData::dimension;

    ParticlesData param;
    auto& particle    = param.particle;
    auto& sourceData  = param.sourceData;
    auto& cellOverlap = param.cellOverlap;
    auto& destData    = param.destData;

    particle.iCell = ConstArray<int, dim>(15);

    sourceData.domainParticles.push_back(particle);

    SAMRAI::tbox::MessageStream particlesWriteStream;

    sourceData.packStream(particlesWriteStream, *cellOverlap);

    SAMRAI::tbox::MessageStream particlesReadStream{particlesWriteStream.getCurrentSize(),
                                                    SAMRAI::tbox::MessageStream::Read,
                                                    particlesWriteStream.getBufferStart()};

    destData.unpackStream(particlesReadStream, *cellOverlap);

    auto expectediCell = ConstArray<int, dim>(-1);

    EXPECT_THAT(destData.patchGhostParticles[0].v, Eq(particle.v));
    EXPECT_THAT(destData.patchGhostParticles[0].iCell, Eq(expectediCell));
    EXPECT_THAT(destData.patchGhostParticles[0].delta, Eq(particle.delta));
    EXPECT_THAT(destData.patchGhostParticles[0].weight, Eq(particle.weight));
    EXPECT_THAT(destData.patchGhostParticles[0].charge, Eq(particle.charge));
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
