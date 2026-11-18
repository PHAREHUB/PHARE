// Fully kinetic: the particles that left the level from a patch (bandLeavers) must reach the
// level ghost particles of every same-level patch whose ghost layer holds them, and nothing
// else must move (docs/superpowers/specs/p4/sawtooth_fix_spec.md, hand-over design).
// Two 2D patches side by side in x make the level: [0,9]x[0,9] (dst) and [10,19]x[0,9] (src),
// periodic in x with period 20.

#include "phare_core.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "simulator/simulator_def.hpp"

#include "amr/data/particles/particles_data.hpp"
#include "amr/data/particles/particles_variable_fill_pattern.hpp"

#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include <set>
#include <array>


using namespace PHARE;

static constexpr std::size_t dim = 2;
using GridLayout_t               = core::PHARE_Types<SimOpts{
                  .dimension = dim, .interp_order = 1, .fully_kinetic = true}>::Hybrid::GridLayout_t;
using ParticleArray_t            = core::ParticleArray<dim>;
using ParticlesData_t            = amr::ParticlesData<ParticleArray_t>;
using Cell                       = std::array<int, dim>;

static_assert(GridLayout_t::options.particle_ghost_width == 3);


SAMRAI::hier::Box box(Cell const& lower, Cell const& upper)
{
    return amr::samrai_box_from(core::Box<int, dim>{lower, upper});
}

std::multiset<Cell> cellsOf(ParticleArray_t const& particles)
{
    std::multiset<Cell> cells;
    for (auto const& particle : particles)
        cells.insert(particle.iCell);
    return cells;
}


struct BandLeaverHandOver : public ::testing::TestWithParam<bool>
{
    SAMRAI::hier::IntVector const ghosts{SAMRAI::tbox::Dimension{dim},
                                         GridLayout_t::options.particle_ghost_width};
    SAMRAI::hier::Box const dstBox = box({0, 0}, {9, 9});
    SAMRAI::hier::Box const srcBox = box({10, 0}, {19, 9});

    ParticlesData_t dst{dstBox, ghosts, "pop"};
    ParticlesData_t src{srcBox, ghosts, "pop"};

    amr::ParticleBandLeaverFillPattern<GridLayout_t> pattern;

    void addTo(ParticleArray_t& particles, Cell const& cell)
    {
        core::Particle<dim> particle;
        particle.weight = 1.;
        particle.iCell  = cell;
        particles.push_back(particle);
    }

    // a periodic shift puts the source on the left of dst
    void handOver(int const xOffset)
    {
        SAMRAI::hier::IntVector offset{SAMRAI::tbox::Dimension{dim}, 0};
        offset[0] = xOffset;
        SAMRAI::hier::Transformation const transformation{offset};

        SAMRAI::pdat::CellGeometry const dstGeometry{dstBox, ghosts};
        SAMRAI::pdat::CellGeometry const srcGeometry{srcBox, ghosts};
        auto fillBox = dstBox;
        fillBox.grow(ghosts);

        auto const overlap = pattern.calculateOverlap(dstGeometry, srcGeometry, dstBox, srcBox,
                                                      fillBox, true, transformation);
        ASSERT_TRUE(dynamic_cast<amr::ParticlesBandOverlap const*>(overlap.get()));

        bool const throughStream = GetParam();
        if (!throughStream)
            dst.copy(src, *overlap);
        else
        {
            SAMRAI::tbox::MessageStream writeStream;
            src.packStream(writeStream, *overlap);
            EXPECT_EQ(writeStream.getCurrentSize(), src.getDataStreamSize(*overlap));
            SAMRAI::tbox::MessageStream readStream{writeStream.getCurrentSize(),
                                                   SAMRAI::tbox::MessageStream::Read,
                                                   writeStream.getBufferStart()};
            dst.unpackStream(readStream, *overlap);
        }
    }

    void expectOnlyLevelGhosts(std::multiset<Cell> const& expected)
    {
        EXPECT_EQ(cellsOf(dst.levelGhostParticles), expected);
        EXPECT_EQ(dst.domainParticles.size(), 0u);
        EXPECT_EQ(dst.patchGhostParticles.size(), 0u);
        EXPECT_EQ(dst.bandLeavers.size(), 0u);
    }
};


TEST_P(BandLeaverHandOver, leaversInDstGhostLayerReachItsLevelGhosts)
{
    addTo(src.bandLeavers, {10, -1}); // below the level, in dst's ghost layer
    addTo(src.bandLeavers, {12, -3}); // corner of dst's ghost layer
    addTo(src.bandLeavers, {13, -1}); // beyond dst's ghost layer
    addTo(src.bandLeavers, {15, 10}); // above the level, far from dst

    // the other arrays of the source never move with this overlap
    addTo(src.domainParticles, {10, 0});
    addTo(src.patchGhostParticles, {9, 0});
    addTo(src.levelGhostParticles, {10, -2});

    handOver(0);

    expectOnlyLevelGhosts({{10, -1}, {12, -3}});
}


TEST_P(BandLeaverHandOver, nothingReachesDstDomain)
{
    addTo(src.bandLeavers, {9, 5}); // not a real leaver: in dst's domain
    handOver(0);
    expectOnlyLevelGhosts({});
}


TEST_P(BandLeaverHandOver, periodicImageLeaversReachDstGhostLayer)
{
    addTo(src.bandLeavers, {19, -1}); // image (-1, -1): dst's ghost layer
    addTo(src.bandLeavers, {17, 10}); // image (-3, 10): dst's ghost layer corner
    addTo(src.bandLeavers, {10, -1}); // image (-10, -1): far from dst

    handOver(-20);

    expectOnlyLevelGhosts({{-1, -1}, {-3, 10}});
}


INSTANTIATE_TEST_SUITE_P(LocalAndStream, BandLeaverHandOver, ::testing::Values(false, true));


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
