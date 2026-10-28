
#include "core/vector.hpp"
#include "core/utilities/kernels.hpp"
#include "core/data/particles/particle.hpp"
#include "core/numerics/pusher/boris/basics.hpp"

#include "gtest/gtest.h"

namespace PHARE::core
{
using Particle_t                        = Particle<1>;
std::array<double, 1> static const half = {1};

auto static make_particle(double const v)
{
    return Particle_t{1., {5}, {.5}, {v, 0, 0}};
}

TEST(BorisMoveError, cpuThrowsWithoutErrorRecord)
{
    auto particle = make_particle(10);
    EXPECT_THROW((boris::advance<AllocatorMode::CPU>(particle, half)), boris::MoveTwoCellException);
}

TEST(BorisMoveError, recordsTheFirstMoveAndKeepsTheParticleInPlace)
{
    boris::MoveError error;
    auto fast   = make_particle(10);
    auto faster = make_particle(20);

    EXPECT_EQ(boris::advance<AllocatorMode::GPU_UNIFIED>(fast, half, &error), fast.iCell());
    EXPECT_EQ(boris::advance<AllocatorMode::CPU>(faster, half, &error), faster.iCell());

    EXPECT_EQ(error.flag, 1);
    EXPECT_DOUBLE_EQ(error.delta, 10.5);
    EXPECT_DOUBLE_EQ(error.vel, 10);
}

TEST(BorisMoveError, noRecordWithinOneCell)
{
    boris::MoveError error;
    auto particle = make_particle(1);
    EXPECT_EQ(boris::advance<AllocatorMode::GPU_UNIFIED>(particle, half, &error)[0], 6);
    EXPECT_EQ(error.flag, 0);
}

#if PHARE_HAVE_MKN_GPU
TEST(BorisMoveError, recordsFromAKernel)
{
    auto errors    = PHARE::Vector<boris::MoveError, AllocatorMode::GPU_UNIFIED>::make(1);
    auto particles = PHARE::Vector<Particle_t, AllocatorMode::GPU_UNIFIED>::make(64);
    for (std::size_t i = 0; i < particles.size(); ++i)
        particles[i] = make_particle(i == 42 ? 10 : 1);

    auto* const ps  = particles.data();
    auto* const err = errors.data();
    auto const dt   = half;
    kernel::launch(particles.size(), [=] _PHARE_ALL_FN_() {
        auto const i  = kernel::idx();
        ps[i].iCell() = boris::advance<AllocatorMode::GPU_UNIFIED>(ps[i], dt, err);
    });

    EXPECT_EQ(errors[0].flag, 1);
    EXPECT_DOUBLE_EQ(errors[0].vel, 10);
    EXPECT_EQ(particles[42].iCell()[0], 5);
    EXPECT_EQ(particles[0].iCell()[0], 6);
}
#endif

} // namespace PHARE::core


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
