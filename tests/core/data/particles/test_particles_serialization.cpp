

#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_serializer.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cmath>


namespace PHARE::core
{
std::uint8_t static constexpr cells = 4;
std::uint8_t static constexpr ppc   = 4;


template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
};


struct AParticleArraySerializationTest : public ::testing::Test
{
};

template<typename Param>
struct ParticleArraySerializationTest : public AParticleArraySerializationTest
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static interp      = 1;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;

    auto constexpr static sim_opts = SimOpts{.dimension = dim, .interp_order = interp,
                                             .layout_mode = layout_mode, .alloc_mode = alloc_mode};

    template<auto lm>
    using Particles_t
        = ParticleArray<ParticleArrayOptions{dim, lm, StorageMode::VECTOR, alloc_mode}>;
    using GridLayout_t     = PHARE_Types<sim_opts>::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;
    using ParticleArray_t  = Particles_t<layout_mode>;

    TestGridLayout_t layout{cells};
};




// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

   TestParam<3, LayoutMode::AoS, AllocatorMode::CPU>
  ,TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>

// PHARE_WITH_THRUST(
//  ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::CPU>
// )

>;
// clang-format on


TYPED_TEST_SUITE(ParticleArraySerializationTest, Permutations_t, );

TYPED_TEST(ParticleArraySerializationTest, test_compare_from_disk)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;

    auto particles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(particles, this->layout.AMRBox(), ppc);

    std::string const file_name = "particles.bin";
    serialize_particles(file_name, particles);

    auto loaded = make_particles<ParticleArray_t>(this->layout);
    deserialize_particles<ParticleArray_t>(file_name, loaded);

    EXPECT_EQ(loaded, particles);
}

TYPED_TEST(ParticleArraySerializationTest, test_deserialize_to_soavx)
{
    using ParticleArray_t  = TestFixture::ParticleArray_t;
    using SoAVXParticles_t = TestFixture::template Particles_t<LayoutMode::SoAVX>;

    auto particles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(particles, this->layout.AMRBox(), ppc);

    std::string const file_name = "particles.bin";
    serialize_particles(file_name, particles);

    auto loaded = make_particles<SoAVXParticles_t>(this->layout);
    deserialize_particles<SoAVXParticles_t, ParticleArray_t>(file_name, loaded);

    EXPECT_TRUE(compare_particles(particles, loaded));
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
