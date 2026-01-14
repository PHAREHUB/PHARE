#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"

#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cmath>


namespace PHARE::core
{
std::size_t static constexpr cells = 3;
std::size_t static constexpr ppc   = 10;


template<std::size_t dim, typename ICell>
PHARE::core::Particle<dim> particle(ICell const& icell)
{
    return {/*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ icell,
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
            /*.v      = */ {{.00001, .00001, .00001}}};
}

template<typename ParticleArray_t, typename Box_t>
void add_particles_in(ParticleArray_t& particles, Box_t const& box)
{
    for (auto const& amr_idx : box)
        for (std::size_t i = 0; i < ppc; ++i)
            particles.emplace_back(particle<ParticleArray_t::dimension>(*amr_idx));
}



template<typename ParticleArray_>
struct ParticleArrayConsistencyTest : public ::testing::Test
{
    auto constexpr static dim    = ParticleArray_::dimension;
    auto constexpr static interp = 1;

    using GridLayout_t = TestGridLayout<typename PHARE_Types<SimOpts<>{dim, interp}>::GridLayout_t>;
    using ParticleArray_t = ParticleArray_;

    GridLayout_t layout{cells};
};



using Permutations_t = testing::Types<ParticleArray<1>>;


TYPED_TEST_SUITE(ParticleArrayConsistencyTest, Permutations_t, );



TYPED_TEST(ParticleArrayConsistencyTest, test_is_consistent_after_swap_copy)
{
    using ParticleArray_t     = TestFixture::ParticleArray_t;
    auto static constexpr dim = ParticleArray_t::dimension;

    auto levelGhostParticles = ParticleArray<dim>{this->layout.AMRBox()};
    add_particles_in(levelGhostParticles, this->layout.AMRBox());

    auto levelGhostParticlesNew = ParticleArray<dim>{this->layout.AMRBox()};
    add_particles_in(levelGhostParticlesNew, this->layout.AMRBox());

    auto levelGhostParticlesOld = ParticleArray<dim>{this->layout.AMRBox()};
    add_particles_in(levelGhostParticlesOld, this->layout.AMRBox());

    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
    levelGhostParticlesNew.clear();
    levelGhostParticles = levelGhostParticlesOld;

    EXPECT_EQ(levelGhostParticlesNew.size(), 0);
    EXPECT_EQ(levelGhostParticlesOld.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(levelGhostParticles.size(), ppc * std::pow(cells, TestFixture::dim));

    EXPECT_TRUE(levelGhostParticlesNew.is_consistent());
    EXPECT_TRUE(levelGhostParticlesOld.is_consistent());
    EXPECT_TRUE(levelGhostParticles.is_consistent());
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
