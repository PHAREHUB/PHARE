
#include <type_traits>


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "data/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/particles/particle_array.h"
#include "data/vecfield/vecfield.h"
#include "hybrid/hybrid_quantities.h"

using namespace PHARE;

class theIons : public ::testing::Test
{
protected:
    Ions<IonPopulation<ParticleArray<1>, VecField<NdArrayVector1D<>, HybridQuantity>>> ions{
        "testIons"};


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




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
