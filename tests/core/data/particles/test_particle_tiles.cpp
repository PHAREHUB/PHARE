#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_tile.hpp"
#include "core/data/tiles/tile_set.hpp"

using namespace PHARE::core;

template<std::size_t dimension>
using ParticleTile_t = ParticleTile<ParticleArray<dimension>>;

template<std::size_t dimension>
using TileSet_t = TileSet<ParticleTile_t<dimension>>;

using DimParticleTiles = testing::Types<ParticleTile_t<1>, ParticleTile_t<2>, ParticleTile_t<3>>;

template<typename ParticleTile>
class ParticleTileTest : public ::testing::Test
{
protected:
    ParticleTileTest() {}
};

TYPED_TEST_SUITE(ParticleTileTest, DimParticleTiles);

TYPED_TEST(ParticleTileTest, constructs)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(50)};
    auto const tile_size = PHARE::core::ConstArray<std::size_t, dim>(4);
    TileSet<TypeParam> particleTileSet{box, tile_size};
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
