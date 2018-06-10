
#include <type_traits>


#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/particle_initializers/fluid_particle_initializer.h"
#include "data/particles/particle_array.h"
#include "utilities/function/function.h"
#include "utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE;


double density(double x)
{
    return 1.;
}



std::array<double, 3> bulkVelocity(double x)
{
    return {{1.0, 0., 0.}};
}



std::array<double, 3> thermalvelocity(double x)
{
    return {{0.2, 0.2, 0.2}};
}




class aFluidParticleInitializer1D : public ::testing::Test
{
private:
    using GridLayoutT    = GridLayout<GridLayoutImplYee<1, 1>>;
    using ParticleArrayT = ParticleArray<1>;

public:
    aFluidParticleInitializer1D()
        : layout{{{0.1}}, {{50}}, Point<double, 1>{0.}}
        , initializer{std::make_unique<FluidParticleInitializer<ParticleArrayT, GridLayoutT>>(
              std::make_unique<ScalarFunction<1>>(density),
              std::make_unique<VectorFunction<1>>(bulkVelocity),
              std::make_unique<VectorFunction<1>>(thermalvelocity), 1., nbrParticlesPerCell)}
    {
        //
    }


    GridLayoutT layout;
    ParticleArrayT particles;
    uint32 nbrParticlesPerCell{1000};
    std::unique_ptr<FluidParticleInitializer<ParticleArrayT, GridLayoutT>> initializer;
};




TEST_F(aFluidParticleInitializer1D, loadsTheCorrectNbrOfParticles)
{
    auto nbrCells             = layout.nbrCells();
    auto expectedNbrParticles = nbrParticlesPerCell * nbrCells[0];
    initializer->loadParticles(particles, layout);
    EXPECT_EQ(expectedNbrParticles, particles.size());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
