
#include <type_traits>


#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_utilities.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

using namespace PHARE::core;
using namespace PHARE::initializer;


class AMaxwellianParticleInitializer1D : public ::testing::Test
{
private:
    using GridLayoutT       = GridLayout<GridLayoutImplYee<1, 1>>;
    using ParticleArrayT    = ParticleArray<1>;
    using InitFunctionArray = std::array<InitFunction<1>, 3>;

public:
    AMaxwellianParticleInitializer1D()
        : layout{{{0.1}}, {{50}}, Point{0.}, Box{Point{50}, Point{99}}}
        , initializer{std::make_unique<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>>(
              density, InitFunctionArray{vx, vy, vz}, InitFunctionArray{vthx, vthy, vthz}, 1.,
              nbrParticlesPerCell)}
    {
        //
    }


    GridLayoutT layout;
    ParticleArrayT particles;
    std::uint32_t nbrParticlesPerCell{10000};
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

        if (!((pos[0] > 0.) and (pos[0] < endDomain)))
            std::cout << "position : " << pos[0] << " not in domain (0," << endDomain << ")\n";
        EXPECT_TRUE(pos[0] > 0. && pos[0] < endDomain);
        i++;
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
