

#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "phare_simulator_options.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/particles/particle_array.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cmath>


namespace PHARE::core
{
std::uint32_t static constexpr cells = 22;
std::size_t static constexpr ppc     = 100;


template<auto _opts>
struct TestParam
{
    auto constexpr static opts        = _opts;
    auto constexpr static dim         = opts.dimension;
    auto constexpr static layout_mode = opts.layout_mode;
    auto constexpr static alloc_mode  = opts.alloc_mode;

    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;
};



template<typename TestParam>
struct ParticleArrayConsistencyTest : public ::testing::Test
{
    auto constexpr static dim = TestParam::dim;
    // auto constexpr static interp = 1;

    using PhareTypes       = PHARE_Types<TestParam::opts>;
    using GridLayout_t     = PhareTypes ::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;
    using ParticleArray_t  = TestParam::ParticleArray_t;



    TestGridLayout_t layout{cells};
};


using Permutations_t = testing::Types<
    TestParam<SimOpts{.dimension = 3, .interp_order = 1, .layout_mode = LayoutMode::AoSTS}>>;

TYPED_TEST_SUITE(ParticleArrayConsistencyTest, Permutations_t, );



TYPED_TEST(ParticleArrayConsistencyTest, test_is_consistent_after_swap_copy)
{
    using ParticleArray_t     = TestFixture::ParticleArray_t;
    auto static constexpr dim = ParticleArray_t::dimension;

    auto levelGhostParticles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticles, this->layout.AMRBox(), ppc);
    EXPECT_NO_THROW(levelGhostParticles.check());

    auto levelGhostParticlesNew = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticlesNew, this->layout.AMRBox(), ppc);
    EXPECT_NO_THROW(levelGhostParticlesNew.check());

    auto levelGhostParticlesOld = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticlesOld, this->layout.AMRBox(), ppc);
    EXPECT_NO_THROW(levelGhostParticlesOld.check());

    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
    EXPECT_NO_THROW(levelGhostParticlesOld.check());
    levelGhostParticlesNew.clear();
    levelGhostParticles = levelGhostParticlesOld;
    EXPECT_NO_THROW(levelGhostParticlesOld.check());

    EXPECT_EQ(levelGhostParticlesNew.size(), 0);
    EXPECT_EQ(levelGhostParticlesOld.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(levelGhostParticles.size(), ppc * std::pow(cells, TestFixture::dim));

    EXPECT_NO_THROW(levelGhostParticlesNew.check());
    EXPECT_NO_THROW(levelGhostParticlesOld.check());
    EXPECT_NO_THROW(levelGhostParticles.check());
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
