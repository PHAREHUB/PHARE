
#include <type_traits>


#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "core/data/particles/particle_array.h"
#include "core/data/particles/particle_utilities.h"
#include "core/utilities/box/box.h"
#include "core/utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;
using namespace PHARE::initializer;

static constexpr size_t PARTICLES_PER_CELL = 1000;

double density(double)
{
    return 1.;
}


double vx(double)
{
    return 1.;
}


double vy(double)
{
    return 1.;
}


double vz(double)
{
    return 1.;
}


double vthx(double)
{
    return 1.;
}


double vthy(double)
{
    return 1.;
}



double vthz(double)
{
    return 1.;
}


template<typename ParticleArray, bool prealloc>
struct AMaxwellianParticleInitializer1D
{
    using GridLayoutT    = GridLayout<GridLayoutImplYee<ParticleArray::dim, 1>>;
    using VectorFunction = std::array<ScalarFunction<1>, 3>;
    using ParticleInit   = MaxwellianParticleInitializer<ParticleArray, GridLayoutT, prealloc>;

    AMaxwellianParticleInitializer1D(size_t alloc = 0)
        : layout{{{0.1}}, {{50}}, Point{0.}, Box{Point{50}, Point{99}}}
        , particles{layout.nbrCells()[0] * alloc}
        , initializer{std::make_unique<ParticleInit>(density, VectorFunction{vx, vy, vz},
                                                     VectorFunction{vthx, vthy, vthz}, 1.,
                                                     nbrParticlesPerCell)}
    {
        if constexpr (!prealloc) // hack sorry
            particles.clear();
    }

    GridLayoutT layout;
    ParticleArray particles;
    uint32 nbrParticlesPerCell{PARTICLES_PER_CELL};
    std::unique_ptr<ParticleInit> initializer;
};

// ParticleArray<1, false> = not contiguous
// AMaxwellianParticleInitializer1D<ParticleArray, false> = not preallocated
using ParticleInitializers
    = testing::Types<AMaxwellianParticleInitializer1D<ParticleArray<1, false>, false>,
                     AMaxwellianParticleInitializer1D<ParticleArray<1, false>, true>,
                     AMaxwellianParticleInitializer1D<ParticleArray<1, true>, false>,
                     AMaxwellianParticleInitializer1D<ParticleArray<1, true>, true>>;

template<typename MaxwellianParticleInitializer>
struct AMaxwellianParticleInitializerTest : public ::testing::Test
{
};

TYPED_TEST_SUITE(AMaxwellianParticleInitializerTest, ParticleInitializers);

TYPED_TEST(AMaxwellianParticleInitializerTest, loadsTheCorrectNbrOfParticles_)
{
    TypeParam init{PARTICLES_PER_CELL};
    auto& layout      = init.layout;
    auto& initializer = init.initializer;
    auto& particles   = init.particles;

    auto expectedNbrParticles = init.nbrParticlesPerCell * layout.nbrCells()[0];
    initializer->loadParticles(particles, layout);
    EXPECT_EQ(expectedNbrParticles, particles.size());
}

TYPED_TEST(AMaxwellianParticleInitializerTest, loadsParticlesInTheDomain_)
{
    TypeParam init{PARTICLES_PER_CELL};
    auto& layout      = init.layout;
    auto& initializer = init.initializer;
    auto& particles   = init.particles;

    initializer->loadParticles(particles, layout);
    for (auto const& particle : particles)
    {
        EXPECT_TRUE(particle.iCell[0] >= 50 && particle.iCell[0] <= 99);
        auto pos       = positionAsPoint(particle, layout);
        auto endDomain = layout.origin()[0] + layout.nbrCells()[0] * layout.meshSize()[0];

        EXPECT_TRUE(pos[0] > 0. && pos[0] < endDomain);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
