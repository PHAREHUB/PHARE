#include <type_traits>


#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_utilities.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include "gtest/gtest.h"

#include "tests/initializer/init_functions.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

using namespace PHARE::core;

namespace PHARE::initializer::test_fn::func_1d
{


class AMaxwellianParticleInitializer1D : public ::testing::Test
{
private:
    using ParticleArrayT    = ParticleArray<1>;
    using InitFunctionArray = std::array<InitFunction<1>, 3>;
    using GridLayoutT       = GridLayout<GridLayoutImplYee<1, 1>>;

public:
    AMaxwellianParticleInitializer1D()
        : layout{{{0.1}}, {{50}}, Point{0.}, Box{Point{50}, Point{99}}}
        , particles{layout.AMRBox()}
        , initializer{std::make_unique<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>>(
              density, InitFunctionArray{vx, vy, vz}, InitFunctionArray{vthx, vthy, vthz}, 1.,
              nbrParticlesPerCell)}
    {
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
    for (auto const& particle : particles)
    {
        EXPECT_TRUE(particle.iCell[0] >= 50 && particle.iCell[0] <= 99);
        auto pos       = positionAsPoint(particle, layout);
        auto endDomain = layout.origin()[0] + layout.nbrCells()[0] * layout.meshSize()[0];

        if (!((pos[0] > 0.) and (pos[0] < endDomain)))
            std::cout << "position : " << pos[0] << " not in domain (0," << endDomain << ")\n";
        EXPECT_TRUE(pos[0] > 0. && pos[0] < endDomain);
    }
}

} // namespace PHARE::initializer::test_fn::func_1d


namespace PHARE::initializer::test_fn::func_2d
{

class AMaxwellianParticleInitializer2D : public ::testing::Test
{
private:
    using ParticleArrayT    = ParticleArray<2>;
    using InitFunctionArray = std::array<InitFunction<2>, 3>;
    using GridLayoutT       = GridLayout<GridLayoutImplYee<2, 1>>;



public:
    AMaxwellianParticleInitializer2D()
        : layout{50}
        , initializer{std::make_unique<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>>(
              density, InitFunctionArray{vx, vy, vz}, InitFunctionArray{vthx, vthy, vthz}, 1.,
              nbrParticlesPerCell)}
    {
    }

    TestGridLayout<GridLayoutT> layout;
    ParticleArrayT particles{layout.AMRBox()};
    std::uint32_t nbrParticlesPerCell{600};
    std::unique_ptr<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>> initializer;
};

TEST_F(AMaxwellianParticleInitializer2D, loadsTheCorrectNbrOfParticles)
{
    // vector push back allocation observations
    // 100 ppc = 262144 - 250000 = 12144 == 12144 * 64 / 1e6 == .7MB overallocated
    // 600 ppc = 2097152 - 1500000 = 597152 * 64 / 1e6 == 38MB overallocated

    auto const expectedNbrParticles = nbrParticlesPerCell * product(layout.AMRBox().shape());
    initializer->loadParticles(particles, layout);
    EXPECT_EQ(expectedNbrParticles, particles.size());
    auto outer_cell_count = std::pow(50, 2) - std::pow(48, 2);
    EXPECT_EQ(particles.capacity(),
              particles.size() + (outer_cell_count * nbrParticlesPerCell * .1));

    // new method
    // 100 ppc = (1511760 - 1500000) * 64 / 1e6 == 0.12544 overallocated
    // 600 ppc = (1511760 - 1500000) * 64 / 1e6 == 0.75264 overallocated
}


} // namespace PHARE::initializer::test_fn::func_2d

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
