//
//
//
//

#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/core/data/electromag/test_electromag_fixtures.hpp"

#include "gtest/gtest.h"

template<typename ParticlesData>
struct PusherTest : public ::testing::Test, public ParticlesData
{
};

using namespace PHARE::core;

template<std::size_t dim>
using AoSGPUParticleArray = ParticleArray<ParticleArrayOptions{
    dim, LayoutMode::AoS, StorageMode::VECTOR, PHARE::AllocatorMode::GPU_UNIFIED}>;

template<std::size_t dim>
using SoAGPUParticleArray = ParticleArray<ParticleArrayOptions{
    dim, LayoutMode::SoA, StorageMode::VECTOR, PHARE::AllocatorMode::GPU_UNIFIED}>;

using ParticleList = testing::Types<AoSParticleArray<3>, AoSGPUParticleArray<3>,
                                    SoAParticleArray<3>, SoAGPUParticleArray<3> //
                                    >;

TYPED_TEST_SUITE(PusherTest, ParticleList, );

TYPED_TEST(PusherTest, doesMoveParticles)
{
    std::size_t constexpr static dim     = TypeParam::dimension;
    std::size_t constexpr static interp  = 1;
    std::size_t constexpr static ppc     = 22;
    std::uint32_t constexpr static cells = 30;
    auto constexpr static no_op          = [](auto& particleRange) { return particleRange; };

    using ParticleArray_t   = TypeParam;
    using PHARE_Types       = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>;
    using GridLayout_t      = TestGridLayout<typename PHARE_Types::GridLayout_t>;
    using Electromag_t      = UsableElectromag<GridLayout_t, ParticleArray_t::alloc_mode>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using ParticleRange     = PHARE::core::IndexRange<ParticleArray_t>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BorisPusher_t = PHARE::core::BorisPusher<dim, ParticleRange, typename Electromag_t::Super,
                                                   Interpolator, BoundaryCondition, GridLayout_t>;

    GridLayout_t layout{cells};
    Electromag_t em{layout};
    ParticleArray_t domainParticles;
    add_particles_in(domainParticles, layout.AMRBox(), ppc);
    shuffle(domainParticles, 1337);
    assert(domainParticles.size() == layout.AMRBox().size() * ppc);

    BorisPusher_t pusher{layout.meshSize(), .001};

    PHARE_LOG_LINE_STR(domainParticles.v(0)[0]);
    PHARE_LOG_LINE_STR(domainParticles.delta(0)[0]);

    Interpolator interpolator;
    auto range = PHARE::core::makeIndexRange(domainParticles);
    pusher.move(
        /*ParticleRange const&*/ range,
        /*ParticleRange&      */ range,
        /*Electromag const&   */ *em,
        /*double mass         */ 1,
        /*Interpolator&*/ interpolator,
        /*GridLayout const&*/ layout, //
        no_op, no_op);

    PHARE_LOG_LINE_STR(domainParticles.v(0)[0]);
    PHARE_LOG_LINE_STR(domainParticles.delta(0)[0]);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
