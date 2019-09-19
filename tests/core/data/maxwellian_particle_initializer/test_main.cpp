
#include <type_traits>


#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data/particles/particle_array.h"
#include "data/particles/particle_utilities.h"
#include "utilities/box/box.h"
#include "utilities/function/function.h"
#include "utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;


double density([[maybe_unused]] double x)
{
    return 1.;
}

std::array<double, 3> bulkVelocity([[maybe_unused]] double x)
{
    return {{1.0, 0., 0.}};
}



std::array<double, 3> thermalvelocity([[maybe_unused]] double x)
{
    return {{0.2, 0.2, 0.2}};
}


class AMaxwellianParticleInitializer1D : public ::testing::Test
{
private:
    using GridLayoutT    = GridLayout<GridLayoutImplYee<1, 1>>;
    using ParticleArrayT = ParticleArray<1>;

public:
    AMaxwellianParticleInitializer1D()
        : layout{{{0.1}}, {{50}}, Point{0.}, Box{Point{50}, Point{99}}}
        , initializer{std::make_unique<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>>(
              density, bulkVelocity, thermalvelocity, 1., nbrParticlesPerCell)}
    {
        //
    }


    GridLayoutT layout;
    ParticleArrayT particles;
    uint32 nbrParticlesPerCell{1000};
    std::unique_ptr<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>> initializer;
};




TEST_F(AMaxwellianParticleInitializer1D, loadsTheCorrectNbrOfParticles)
{
    auto nbrCells             = layout.nbrCells();
    auto expectedNbrParticles = nbrParticlesPerCell * nbrCells[0];
    initializer->loadParticles(particles, layout);
    EXPECT_EQ(expectedNbrParticles, particles.size());
}




TEST_F(AMaxwellianParticleInitializer1D, loadsParticlesInTheDomain)
{
    initializer->loadParticles(particles, layout);
    auto i = 0u;
    for (auto const& particle : particles)
    {
        EXPECT_TRUE(particle.iCell[0] >= 50 && particle.iCell[0] <= 99);
        auto pos       = positionAsPoint(particle, layout);
        auto endDomain = layout.origin()[0] + layout.nbrCells()[0] * layout.meshSize()[0];

        EXPECT_TRUE(pos[0] > 0. && pos[0] < endDomain);
        i++;
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
