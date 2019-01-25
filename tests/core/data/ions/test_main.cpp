
#include <type_traits>


#include "data/ions/ion_initializer.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "hybrid/hybrid_quantities.h"


#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/particle_initializers/fluid_particle_initializer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using FluidParticleInitializer1D         = FluidParticleInitializer<ParticleArray<1>, GridYee1D>;



double density(double x)
{
    return x * x + 2.;
}

std::array<double, 3> bulkVelocity(double x)
{
    return std::array<double, 3>{{1.0, 0.0, 0.0}};
}


std::array<double, 3> thermalVelocity(double x)
{
    return std::array<double, 3>{{0.5, 0.0, 0.0}};
}




class theIons : public ::testing::Test
{
protected:
    using VecField1D = VecField<NdArrayVector1D<>, HybridQuantity>;

    using IonPopulation1D = IonPopulation<ParticleArray<1>, VecField1D, GridYee1D>;
    Ions<IonPopulation1D, GridYee1D> ions;

    IonsInitializer<ParticleArray<1>, GridYee1D> createInitializer()
    {
        IonsInitializer<ParticleArray<1>, GridYee1D> initializer;

        initializer.masses.push_back(1.);
        initializer.names.push_back("protons");
        initializer.nbrPopulations = 1;
        initializer.name           = "TestIons";
        initializer.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
            std::make_unique<ScalarFunction<1>>(density),
            std::make_unique<VectorFunction<1>>(bulkVelocity),
            std::make_unique<VectorFunction<1>>(thermalVelocity), -1., 10));

        return initializer;
    }


    theIons()
        : ions{createInitializer()}
    {
    }

public:
    ~theIons();
};

theIons::~theIons() {}


TEST_F(theIons, areAContainerOfIonPopulations)
{
    //
    for (auto& pop : ions)
    {
        (void)pop;
    }
}




TEST_F(theIons, areNotUsableUponConstruction)
{
    EXPECT_FALSE(ions.isUsable());
}




TEST_F(theIons, areSettableUponConstruction)
{
    EXPECT_TRUE(ions.isSettable());
}




TEST_F(theIons, throwIfAccessingDensityWhileNotUsable)
{
    EXPECT_ANY_THROW(ions.density());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
