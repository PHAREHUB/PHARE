
#include "core/def/phare_mpi.hpp"


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


template<typename dimType>
struct AParticlesDataND : public testing::Test
{
    static constexpr auto dim = dimType{}();

    SAMRAI::tbox::Dimension dimension{dim};
    SAMRAI::hier::BlockId blockId{0};

    SAMRAI::hier::Box sourceDomain{SAMRAI::hier::Index{dimension, 6},
                                   SAMRAI::hier::Index{dimension, 11}, blockId};

    SAMRAI::hier::Box destDomain{SAMRAI::hier::Index{dimension, 1},
                                 SAMRAI::hier::Index{dimension, 5}, blockId};

    SAMRAI::hier::IntVector ghost{SAMRAI::hier::IntVector::getOne(dimension)};

    ParticlesData<ParticleArray<dim>> destData{destDomain, ghost, "name"};
    ParticlesData<ParticleArray<dim>> sourceData{sourceDomain, ghost, "name"};
    typename ParticleArray<dim>::Particle_t particle;

    AParticlesDataND()
    {
        particle.weight = 1.0;
        particle.charge = 1.0;
        particle.v      = {1.0, 1.0, 1.0};
    }
};



using WithAllDim = testing::Types<DimConst<1>, DimConst<2>, DimConst<3>>;

TYPED_TEST_SUITE(AParticlesDataND, WithAllDim);


TYPED_TEST(AParticlesDataND, copy_test)
{
    // particle is in the ghost of the source patchdata
    static constexpr auto dim = TypeParam{}();

    this->particle.iCell = ConstArray<int, dim>(5);

    this->sourceData.patchGhostParticles.push_back(this->particle);
    this->destData.copy(this->sourceData);

    ASSERT_THAT(this->destData.domainParticles.size(), Eq(1));
    ASSERT_THAT(this->destData.patchGhostParticles.size(), Eq(0));
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
