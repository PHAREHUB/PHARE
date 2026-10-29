
#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE;
using namespace PHARE::core;

template<typename Simulator>
struct ParticleListTest : public ::testing::Test
{
};

using ParticleList = testing::Types<Particle<1>, Particle<2>, Particle<3>>;

TYPED_TEST_SUITE(ParticleListTest, ParticleList);

TYPED_TEST(ParticleListTest, SoAandAoSInterop)
{
    using Particle             = TypeParam;
    constexpr auto dim         = Particle::dimension;
    constexpr std::size_t size = 10;
    constexpr Box<int, dim> domain{ConstArray<int, dim>(0), ConstArray<int, dim>(size - 1)};

    using Hybrid        = PHARE_Types<SimOpts{dim, 1}>::Hybrid;
    using ParticleAoS_t = Hybrid::ParticleAoS_t;
    using ParticleSoA_t = Hybrid::ParticleSoA_t;

    ParticleSoA_t soa{size};
    for (std::size_t i = 0; i < size; i++)
    {
        soa.weight(i) = 1 + i;
        soa.iCell(i)  = ConstArray<int, dim>(i);
        soa.delta(i)  = ConstArray<double, dim>(i + 1);
        soa.v(i)      = ConstArray<double, 3>(soa.weight(i) + 2);
    }
    EXPECT_EQ(soa.size(), size);

    for (std::size_t i = 0; i < size; i++)
    {
        EXPECT_EQ(soa.weight()[i], i + 1);
        EXPECT_EQ(soa.weight(i), i + 1);
    }

    ParticleAoS_t aos{domain};
    for (std::size_t i = 0; i < size; i++)
    {
        aos.emplace_back(soa.copy(i));
        EXPECT_EQ(soa.copy(i), aos.back());
    }
    EXPECT_EQ(aos.size(), size);
    EXPECT_EQ(soa.size(), aos.size());

    ParticleSoA_t soaFromAoS{aos.size()};
    ParticlePacker<ParticleAoS_t>{aos}.pack(soaFromAoS);

    for (std::size_t i = 0; i < size; i++)
        EXPECT_EQ(soaFromAoS.copy(i), aos[i]);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
