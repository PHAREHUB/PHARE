

#include "core/utilities/types.hpp"
#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_packer.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

template<typename Simulator>
struct ParticleListTest : public ::testing::Test
{
};

// using ParticleList = testing::Types<Particle<1>, Particle<2>, Particle<3>>;

// TYPED_TEST_SUITE(ParticleListTest, ParticleList);

// TYPED_TEST(ParticleListTest, SoAandAoSInterop)
// {
//     using Particle             = TypeParam;
//     constexpr auto dim         = Particle::dimension;
//     constexpr std::size_t size = 10;

//     ParticleArray<dim, true> contiguous{size};
//     for (std::size_t i = 0; i < size; i++)
//     {
//         contiguous.weight[i] = 1 + i;
//         contiguous.charge[i] = 1 + i;
//         contiguous.iCell[i]  = ConstArray<int, dim>(i);
//         contiguous.delta[i]  = ConstArray<double, dim>(i + 1);
//         contiguous.v[i]      = ConstArray<double, 3>(contiguous.weight[i] + 2);
//         EXPECT_EQ(std::copy(contiguous[i]), contiguous[i]);
//     }
//     EXPECT_EQ(contiguous.size(), size);

//     for (std::size_t i = 0; i < size; i++)
//     {
//         EXPECT_EQ(contiguous.weight[i], i + 1); // fastest
//         EXPECT_EQ(contiguous[i].weight, i + 1);
//         EXPECT_EQ(contiguous[i], std::copy(contiguous[i]));
//     }

//     ParticleArray<dim> particleArray;
//     for (auto const& view : contiguous)
//     {
//         auto i = particleArray.size();
//         particleArray.emplace_back(std::copy(view));
//         EXPECT_EQ(contiguous[i], particleArray.back());
//     }
//     EXPECT_EQ(particleArray.size(), size);
//     EXPECT_EQ(contiguous.size(), particleArray.size());

//     ParticleArray<dim, true> AoSFromSoA{particleArray.size()};
//     ParticlePacker<dim>{particleArray}.pack(AoSFromSoA);

//     std::size_t i = 0;
//     for (auto const& particle : AoSFromSoA)
//         EXPECT_EQ(particle, particleArray[i++]);
// }
