
#include "phare_core.hpp"
#include "core/utilities/types.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_comparator.hpp"
#include "core/data/particles/particle_array_converter.hpp"

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


// 2^dim "corner" neighbours only (every coordinate moves by +/-1), ordered so that
// treating -1 as bit 0 and +1 as bit 1 counts up from (-1,-1,-1) to (1,1,1).
template<std::size_t dim>
auto corner_offsets()
{
    std::array<std::array<int, dim>, (std::size_t{1} << dim)> offsets{};
    for (std::size_t idx = 0; idx < offsets.size(); ++idx)
        for (std::size_t d = 0; d < dim; ++d)
            offsets[idx][d] = ((idx >> (dim - 1 - d)) & 1) ? 1 : -1;
    return offsets;
}

template<typename ICell>
auto add_icell(ICell const& a, ICell const& b)
{
    ICell out;
    for (std::size_t d = 0; d < a.size(); ++d)
        out[d] = a[d] + b[d];
    return out;
}


// Function templates can't be partially specialised (fixing layout_mode while leaving
// Particles deduced), so the per-layout dispatch lives on a class template instead,
// which can be. One specialisation per layout_mode, each still generic over Particles.
template<auto layout_mode>
struct MoveParticles;

template<>
struct MoveParticles<LayoutMode::AoS>
{
    // no registration needed, flat array has no cell/tile structure to maintain
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        std::size_t counter = 0;
        for (auto& p : particles)
        {
            p.iCell() = add_icell(p.iCell(), offsets[counter % offsets.size()]);
            ++counter;
        }
    }
};

template<>
struct MoveParticles<LayoutMode::AoSMapped>
{
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        std::size_t counter = 0;
        for (std::size_t idx = 0; idx < particles.size(); ++idx)
        {
            auto const newcell
                = add_icell(particles[idx].iCell(), offsets[counter % offsets.size()]);
            ++counter;
            particles.change_icell(newcell, idx);
        }
    }
};

template<>
struct MoveParticles<LayoutMode::AoSPC>
{
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        for (auto const& bix : particles.local_box())
        {
            auto& cell_particles = particles(bix);
            auto const n         = cell_particles.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto& p            = cell_particles[i];
                auto const oldcell = p.iCell();
                auto const newcell = add_icell(oldcell, offsets[counter % offsets.size()]);
                ++counter;
                p.iCell() = newcell;
                ParticleTracker<dim> const pt{oldcell};
                particles.template move_check<ParticleType::Domain>(pt, i, p);
            }
        }
        sync_moved_pc<ParticleType::Domain>(particles);
    }
};

template<>
struct MoveParticles<LayoutMode::AoSTS>
{
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        for (auto& tile : particles())
        {
            auto& tile_particles = tile();
            auto const n         = tile_particles.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto& p              = tile_particles[i];
                auto const oldcell   = p.iCell();
                auto const tile_cell = particles.local_tile_cell(oldcell);
                auto const newcell   = add_icell(oldcell, offsets[counter % offsets.size()]);
                ++counter;
                p.iCell() = newcell;
                // production contract (see MultiBoris): only particles staying in the patch
                // box register; patch-leavers stay put in their tile with a ghost iCell
                if (isIn(newcell, particles.box()))
                {
                    TiledParticleTracker<dim> const pt{oldcell, tile_cell};
                    particles.template move_check<ParticleType::Domain>(pt, i, p);
                }
            }
        }
        sync_ts<ParticleType::Domain>(particles);
    }
};

template<>
struct MoveParticles<LayoutMode::AoSPCTS>
{
    // registration goes through the span (move_check lives on PCTileSetSpan); the view's
    // gap arrays alias the vector's, so sync_pc_ts on the vector sees every registration
    template<typename Particles, typename Offsets>
    static void apply(Particles& particles, Offsets const& offsets)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        auto view           = particles.view();
        for (auto& tile : view())
        {
            auto& tile_particles = tile();
            for (auto const& bix : tile_particles.local_box())
            {
                auto& cell_particles = tile_particles(bix);
                auto const n         = cell_particles.size();
                for (std::size_t i = 0; i < n; ++i)
                {
                    auto& p              = cell_particles[i];
                    auto const oldcell   = p.iCell();
                    auto const tile_cell = view.local_tile_cell(oldcell);
                    auto const newcell   = add_icell(oldcell, offsets[counter % offsets.size()]);
                    ++counter;
                    p.iCell() = newcell;
                    TiledParticleTracker<dim> const pt{oldcell, tile_cell};
                    view.template move_check<ParticleType::Domain>(pt, i, p);
                }
            }
        }
        sync_pc_ts<ParticleType::Domain>(particles);
    }
};


template<typename Particles>
void move_particles(Particles& particles)
{
    MoveParticles<Particles::layout_mode>::apply(particles, corner_offsets<Particles::dimension>());
}


// level ghost arrays: particles live in the ghost layer; movers re-bucket anywhere
// inside the ghost box (including into the domain), ghost-box leavers are deleted at
// sync. only per-cell layouts register level ghost moves — other layouts skip the test
template<auto layout_mode>
struct MoveLevelGhostParticles;

template<>
struct MoveLevelGhostParticles<LayoutMode::AoS>
{
    // no registration: move directly, then apply the deletion contract by hand
    template<typename Particles, typename Offsets, typename Box_t>
    static void apply(Particles& particles, Offsets const& offsets, Box_t const& ghost_box)
    {
        std::size_t counter = 0;
        for (auto& p : particles)
        {
            p.iCell() = add_icell(p.iCell(), offsets[counter % offsets.size()]);
            ++counter;
        }
        std::erase_if(particles.vector(),
                      [&](auto const& p) { return not isIn(p.iCell(), ghost_box); });
    }
};

template<>
struct MoveLevelGhostParticles<LayoutMode::AoSPC>
{
    template<typename Particles, typename Offsets, typename Box_t>
    static void apply(Particles& particles, Offsets const& offsets, Box_t const&)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        for (auto const& bix : particles.local_box())
        {
            auto& cell_particles = particles(bix);
            auto const n         = cell_particles.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto& p            = cell_particles[i];
                auto const oldcell = p.iCell();
                p.iCell()          = add_icell(oldcell, offsets[counter % offsets.size()]);
                ++counter;
                ParticleTracker<dim> const pt{oldcell};
                particles.template move_check<ParticleType::LevelGhost>(pt, i, p);
            }
        }
        sync_moved_pc<ParticleType::LevelGhost>(particles);
    }
};

template<>
struct MoveLevelGhostParticles<LayoutMode::AoSPCTS>
{
    // registration goes through the span, as for the domain mover above
    template<typename Particles, typename Offsets, typename Box_t>
    static void apply(Particles& particles, Offsets const& offsets, Box_t const&)
    {
        auto constexpr dim  = Particles::dimension;
        std::size_t counter = 0;
        auto view           = particles.view();
        for (auto& tile : view())
        {
            auto& tile_particles = tile();
            for (auto const& bix : tile_particles.local_box())
            {
                auto& cell_particles = tile_particles(bix);
                auto const n         = cell_particles.size();
                for (std::size_t i = 0; i < n; ++i)
                {
                    auto& p              = cell_particles[i];
                    auto const oldcell   = p.iCell();
                    auto const tile_cell = view.local_tile_cell(oldcell);
                    p.iCell() = add_icell(oldcell, offsets[counter % offsets.size()]);
                    ++counter;
                    auto const pt
                        = make_particle_tracker<LayoutMode::AoSPCTS, ParticleType::LevelGhost,
                                                dim>(oldcell, tile_cell);
                    view.template move_check<ParticleType::LevelGhost>(pt, i, p);
                }
            }
        }
        sync_pc_ts<ParticleType::LevelGhost>(particles);
    }
};


// every per-cell bucket must only hold particles whose iCell maps to that bucket
template<typename Particles>
void check_cell_buckets(Particles const& particles)
{
    using enum LayoutMode;
    auto constexpr static dim = Particles::dimension;

    if constexpr (any_in(Particles::layout_mode, AoSPC, AoSPCTS))
    {
        auto const check = [](auto const& cps, std::array<std::uint32_t, dim> const& cell) {
            for (auto const& p : cps(cell))
                EXPECT_TRUE(array_equals(cps.local_cell(p.iCell()), cell))
                    << "particle in wrong cell bucket: " << Point{p.iCell()};
        };

        if constexpr (Particles::layout_mode == AoSPC)
        {
            for (auto const& bix : particles.local_box())
                check(particles, bix);
        }
        else
        {
            for (auto const& tile : particles())
                for (auto const& bix : tile().local_box())
                    check(tile(), bix);
        }
    }
}


// only particles outside the patch domain box may sit outside their tile's box —
// a particle inside the patch domain must never be stored in a tile's ghost cells
template<typename Particles>
void check_tile_ownership(Particles const& particles)
{
    using enum LayoutMode;

    if constexpr (any_in(Particles::layout_mode, AoSTS, AoSPCTS))
    {
        auto const check = [&](auto const& tile, auto const& p) {
            if (isIn(p.iCell(), particles.box()))
                EXPECT_TRUE(isIn(p.iCell(), tile)) << "domain particle in tile ghost cells: "
                                                   << Point{p.iCell()} << " tile: " << tile;
        };

        if constexpr (any_in(Particles::layout_mode, AoSTS))
        {
            for (auto const& tile : particles())
                for (auto const& p : tile())
                    check(tile, p);
        }
        else
        {
            for (auto const& tile : particles())
            {
                auto const& cps = tile();
                for (auto const& bix : cps.local_box())
                    for (auto const& p : cps(bix))
                        check(tile, p);
            }
        }
    }
}



TYPED_TEST(ParticleArrayConstructionTest, test_move_sync_works)
{
    using ParticleArray_t    = TestFixture::ParticleArray_t;
    using AoSParticleArray_t = AoSParticleArray<TestFixture::dim>;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    auto particles = make_particles<ParticleArray_t>(this->layout);
    add_particles_in(particles, this->layout.AMRBox(), ppc);

    // independent AoS ground truth, same initial particles, moved directly (no
    // registration needed for AoS) so it never depends on the layout under test
    auto reference = convert_particles<AoSParticleArray_t>(particles, this->layout);

    move_particles(particles);
    move_particles(reference);

    check_tile_ownership(particles);
    check_cell_buckets(particles);

    // particles have moved out of AMRBox into the first ghost layer; sorting must use the
    // grown box or distinct out-of-box cells alias to the same flat index and, with all
    // deltas identical, the two arrays end up in different orders
    auto const sort_box = grow(this->layout.AMRBox(), 1);
    auto converted      = convert_particles<AoSParticleArray_t>(particles, this->layout);
    sort_particles(converted, sort_box);
    sort_particles(reference, sort_box);

    auto const report = compare_particles(reference, converted);
    EXPECT_TRUE(report) << report.why();
}


// level ghost arrays fill the ghost layer; after a move they may re-bucket anywhere in
// the ghost box (including into the domain) and ghost-box leavers must be deleted
TYPED_TEST(ParticleArrayConstructionTest, test_level_ghost_move_sync_works)
{
    using ParticleArray_t    = TestFixture::ParticleArray_t;
    using AoSParticleArray_t = AoSParticleArray<TestFixture::dim>;
    using enum LayoutMode;

    auto constexpr static dim = TestFixture::dim;

    PHARE_LOG_LINE_SS(ParticleArray_t::type_id);

    if constexpr (ParticleArray_t::alloc_mode != AllocatorMode::CPU
                  or not any_in(ParticleArray_t::layout_mode, AoS, AoSPC, AoSPCTS))
        GTEST_SKIP() << "level ghost move_check unsupported for this layout";
    else
    {
        auto constexpr static ghosts = TestFixture::GridLayout_t::nbrParticleGhosts();
        auto const ghost_box         = grow(this->layout.AMRBox(), ghosts);

        auto particles = make_particles<ParticleArray_t>(this->layout);
        add_ghost_particles(particles, this->layout.AMRBox(), ppc, ghosts);

        auto reference = convert_particles<AoSParticleArray_t>(particles, this->layout);

        auto const offsets = corner_offsets<dim>();
        MoveLevelGhostParticles<ParticleArray_t::layout_mode>::apply(particles, offsets,
                                                                     ghost_box);
        MoveLevelGhostParticles<LayoutMode::AoS>::apply(reference, offsets, ghost_box);

        check_tile_ownership(particles);
        check_cell_buckets(particles);

        auto converted = convert_particles<AoSParticleArray_t>(particles, this->layout);
        sort_particles(converted, ghost_box);
        sort_particles(reference, ghost_box);

        EXPECT_EQ(reference.size(), converted.size());
        auto const report = compare_particles(reference, converted);
        EXPECT_TRUE(report) << report.why();
    }
}


} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
