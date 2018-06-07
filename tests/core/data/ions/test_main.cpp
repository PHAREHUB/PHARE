
#include <type_traits>


#include "data/ions/ion_initializer.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "hybrid/hybrid_quantities.h"



#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE;

struct GridLayoutMock
{
};

class theIons : public ::testing::Test
{
protected:
    Ions<IonPopulation<ParticleArray<1>, VecField<NdArrayVector1D<>, HybridQuantity>>,
         GridLayoutMock>
        ions;
    IonsInitializer<ParticleArray<1>, GridLayoutMock> createInitializer()
    {
        IonsInitializer<ParticleArray<1>, GridLayoutMock> initializer;

        initializer.masses.push_back(1.);
        initializer.names.push_back("protons");
        initializer.nbrPopulations = 1;
        initializer.name           = "TestIons";

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
        //
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



TEST_F(theIons, throwIfAccessingVelocityWhileNotUsable)
{
    EXPECT_ANY_THROW(ions.velocity());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
