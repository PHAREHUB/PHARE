
#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_comparator.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include "gtest/gtest.h"
#include <cmath>


namespace PHARE::core
{
auto static const cells = get_env_as("PHARE_CELLS", std::uint32_t{4});
auto static const ppc   = get_env_as("PHARE_PPC", std::size_t{10});


template<std::size_t _dim, auto lm, auto am>
struct TestParam
{
    static_assert(all_are<LayoutMode>(lm));
    static_assert(all_are<AllocatorMode>(am));

    auto constexpr static dim         = _dim;
    auto constexpr static layout_mode = lm;
    auto constexpr static alloc_mode  = am;
};


template<typename Param>
struct ParticleArrayConstructionTest : public ::testing::Test
{
    auto constexpr static dim         = Param::dim;
    auto constexpr static layout_mode = Param::layout_mode;
    auto constexpr static alloc_mode  = Param::alloc_mode;
    auto constexpr static sim_opts
        = SimOpts{.dimension = dim, .layout_mode = layout_mode, .alloc_mode = alloc_mode};

    using GridLayout_t = TestGridLayout<typename PHARE_Types<sim_opts>::GridLayout_t>;
    using ParticleArray_t
        = ParticleArray<ParticleArrayOptions{dim, layout_mode, StorageMode::VECTOR, alloc_mode}>;

    GridLayout_t layout{cells};

    ParticleArray_t setup_particles() const // test movable
    {
        auto ps = make_particles<ParticleArray_t>(layout);
        add_particles_in(ps, layout.AMRBox(), ppc);
        return ps;
    }
};



template<typename ParticleArrayConstructionTest_t>
auto run(ParticleArrayConstructionTest_t& self)
{
    // abort_if(self.periodic_neighbours_for(13).size());
}

// clang-format off
using Permutations_t = testing::Types< // ! notice commas !

    TestParam<3, LayoutMode::AoS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSMapped, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSPC, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::CPU>
   ,TestParam<3, LayoutMode::AoSPCTS, AllocatorMode::CPU>

PHARE_WITH_THRUST(
    ,TestParam<3, LayoutMode::SoA, AllocatorMode::CPU>
    ,TestParam<3, LayoutMode::SoAPC, AllocatorMode::CPU>
)

PHARE_WITH_GPU(
    ,TestParam<3, LayoutMode::AoS, AllocatorMode::GPU_UNIFIED>
    ,TestParam<3, LayoutMode::AoSTS, AllocatorMode::GPU_UNIFIED>
)

>;
// clang-format on



TYPED_TEST_SUITE(ParticleArrayConstructionTest, Permutations_t, );



TYPED_TEST(ParticleArrayConstructionTest, test_move_swap_etc)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto levelGhostParticles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticles, this->layout.AMRBox(), ppc);

    auto levelGhostParticlesNew = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticlesNew, this->layout.AMRBox(), ppc);

    auto levelGhostParticlesOld = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(levelGhostParticlesOld, this->layout.AMRBox(), ppc);

    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
    levelGhostParticlesNew.clear();
    levelGhostParticles.clear();
    levelGhostParticles = levelGhostParticlesOld;

    EXPECT_EQ(levelGhostParticlesNew.size(), 0);
    EXPECT_EQ(levelGhostParticlesOld.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(levelGhostParticles.size(), ppc * std::pow(cells, TestFixture::dim));

    // size matching isn't enough: assignment (unlike copy construction) must also
    // leave levelGhostParticles' content/views correctly pointing at its own copied
    // data, not stale/aliased to levelGhostParticlesOld's storage
    per_particle(levelGhostParticles, [&](auto const& p) {
        EXPECT_NE(p.v()[0], 0);
        EXPECT_NE(p.v()[1], 0);
        EXPECT_NE(p.v()[2], 0);
    });

    auto const eq = compare_particles(levelGhostParticles, levelGhostParticlesOld);
    if (!eq)
    {
        PHARE_LOG_LINE_SS(eq.why());
    }
    EXPECT_TRUE(eq);

    // the assignment's real risk is a stale view: for VECTOR storage, constructing a
    // span-mode view reads through particles_views_ (rebuilt by sync()/reset_views(),
    // not by operator=), so check the view itself, not just the vector-mode array
    auto view = levelGhostParticles.view();

    per_particle(view, [&](auto const& p) {
        EXPECT_NE(p.v()[0], 0);
        EXPECT_NE(p.v()[1], 0);
        EXPECT_NE(p.v()[2], 0);
    });

    auto const view_eq = compare_particles(view, levelGhostParticlesOld);
    if (!view_eq)
    {
        PHARE_LOG_LINE_SS(view_eq.why());
    }
    EXPECT_TRUE(view_eq);
}



TYPED_TEST(ParticleArrayConstructionTest, test_copy_move)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto particles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(particles, this->layout.AMRBox(), ppc);

    auto const copy = particles;
    auto const move = std::move(particles);

    EXPECT_EQ(move.size(), ppc * std::pow(cells, TestFixture::dim));
    EXPECT_EQ(move.size(), copy.size());

    auto const eq = compare_particles(copy, move);
    if (!eq)
    {
        PHARE_LOG_LINE_SS(eq.why());
    }
    EXPECT_TRUE(eq);
}



TYPED_TEST(ParticleArrayConstructionTest, test_move_on_create)
{
    using ParticleArray_t = TestFixture::ParticleArray_t;

    // static_assert(std::is_trivially_move_assignable_v<ParticleArray_t>);
    // static_assert(std::is_trivially_move_constructible_v<ParticleArray_t>);

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    std::vector<ParticleArray_t> vecs;
    vecs.emplace_back(this->setup_particles());
    ParticleArrayService::sync<0, ParticleType::Domain>(vecs.back());
    check_particles_views(vecs.back());

    // EXPECT_EQ(move.size(), ppc * std::pow(cells, TestFixture::dim));
    // EXPECT_EQ(copy.size(), ppc * std::pow(cells, TestFixture::dim));
    // EXPECT_TRUE(compare_particles(copy, move));
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
