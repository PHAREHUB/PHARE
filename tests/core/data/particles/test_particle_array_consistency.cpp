#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <tuple>
#include <vector>
#include <algorithm>


namespace PHARE::core
{
std::uint32_t static constexpr cells = 14;
std::size_t static constexpr ppc     = 100;


template<std::size_t dim, typename ICell>
PHARE::core::Particle<dim> particle(ICell const& icell)
{
    return {/*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ icell,
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
            /*.v      = */ {{.00001, .00001, .00001}}};
}



template<typename ParticleArray_>
struct ParticleArrayConsistencyTest : public ::testing::Test
{
    auto constexpr static dim    = ParticleArray_::dimension;
    auto constexpr static interp = 1;

    using GridLayout_t
        = TestGridLayout<typename PHARE_Types<SimOpts{dim, interp}>::Hybrid::GridLayout_t>;
    using ParticleArray_t = ParticleArray_;

    GridLayout_t layout{cells};
};


using Permutations_t = testing::Types<ParticleArray<ParticleArrayOptions{1}>,
                                      ParticleArray<ParticleArrayOptions{3, LayoutMode::AoSPCTS}>>;


TYPED_TEST_SUITE(ParticleArrayConsistencyTest, Permutations_t, );


TYPED_TEST(ParticleArrayConsistencyTest, test_is_consistent_after_swap_copy)
{
    using ParticleArray_t     = TestFixture::ParticleArray_t;
    auto static constexpr dim = ParticleArray_t::dimension;

    auto levelGhostParticles = make_particles<ParticleArray_t>(this->layout);
    add_particles(levelGhostParticles, this->layout.AMRBox(), ppc);

    auto levelGhostParticlesNew = make_particles<ParticleArray_t>(this->layout);
    add_particles(levelGhostParticlesNew, this->layout.AMRBox(), ppc);

    auto levelGhostParticlesOld = make_particles<ParticleArray_t>(this->layout);
    add_particles(levelGhostParticlesOld, this->layout.AMRBox(), ppc);

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


// mirrors HybridHybridMessengerStrategy::lastStep() (swap New<->Old, clear New, copy Old
// into the pushable buffer) but with distinguishable per-particle delta/v (rather than the
// fixed constants add_particles() uses) so a value taken from the wrong particle, or
// corrupted in transit, is actually detectable - not just silently identical to its sibling.
TYPED_TEST(ParticleArrayConsistencyTest, test_swap_copy_preserves_values)
{
    using ParticleArray_t     = TestFixture::ParticleArray_t;
    auto static constexpr dim = ParticleArray_t::dimension;

    auto levelGhostParticles    = make_particles<ParticleArray_t>(this->layout);
    auto levelGhostParticlesNew = make_particles<ParticleArray_t>(this->layout);
    auto levelGhostParticlesOld = make_particles<ParticleArray_t>(this->layout);

    // per_particle() (particle_array.hpp:230) visits tile-then-cell for tiled layouts
    // instead of the array's own top-level begin()/end(), which tiled layouts don't support.
    auto disperse = [](auto& arr, int seed) {
        auto gen = rando(seed);
        ParticleDeltaDistribution<double> deltaDistrib;
        std::uniform_real_distribution<double> veloc{-6, 6};
        per_particle(arr, [&](auto& p) {
            p.delta() = core::ConstArrayFrom<dim>([&] { return deltaDistrib(gen); });
            p.v()     = core::ConstArrayFrom<3>([&] { return veloc(gen); });
        });
    };

    add_particles(levelGhostParticlesOld, this->layout.AMRBox(), ppc);
    disperse(levelGhostParticlesOld, 1337);

    add_particles(levelGhostParticlesNew, this->layout.AMRBox(), ppc);
    disperse(levelGhostParticlesNew, 42);

    auto const snapshot = [](auto const& arr) {
        std::vector<Particle<dim>> out;
        per_particle(arr, [&](auto const& p) { out.push_back(p); });
        std::sort(out.begin(), out.end(), [](auto const& a, auto const& b) {
            return std::tie(a.iCell(), a.delta()) < std::tie(b.iCell(), b.delta());
        });
        return out;
    };

    // levelGhostParticlesNew's content is what levelGhostParticles must equal afterwards
    auto const expected = snapshot(levelGhostParticlesNew);

    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
    levelGhostParticlesNew.clear();
    levelGhostParticles = levelGhostParticlesOld;

    auto const actual = snapshot(levelGhostParticles);

    ASSERT_EQ(expected.size(), actual.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        EXPECT_EQ(expected[i].iCell(), actual[i].iCell());
        EXPECT_EQ(expected[i].delta(), actual[i].delta());
        EXPECT_EQ(expected[i].v(), actual[i].v());
    }
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
