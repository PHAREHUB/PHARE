#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "core/data/clusters/clusters.hpp"
#include "core/data/particles/particle_cluster.hpp"
#include "core/data/particles/particle_array.hpp"


using namespace PHARE::core;

template<std::size_t dimension>
using ParticleCluster_t = ParticleCluster<ParticleArray<dimension>>;

template<std::size_t dimension>
using ClusterSet_t = ClusterSet<ParticleCluster_t<dimension>>;

using DimParticleClusters
    = testing::Types<ParticleCluster_t<1>, ParticleCluster_t<2>, ParticleCluster_t<3>>;


template<typename ParticleCluster>
class ParticleClusterTest : public ::testing::Test
{
protected:
    ParticleClusterTest() {}
};

TYPED_TEST_SUITE(ParticleClusterTest, DimParticleClusters);

TYPED_TEST(ParticleClusterTest, constructs)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(50)};
    auto const cluster_size = PHARE::core::ConstArray<int, dim>(4);
    ClusterSet<TypeParam> particleClusterSet{box, cluster_size};
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
