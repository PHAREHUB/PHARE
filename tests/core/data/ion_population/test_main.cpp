
#include <type_traits>




#include "data/ions/ion_population/ion_population.h"
#include "data/particles/particle_array.h"
#include "hybrid/hybrid_quantities.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;




struct DummyField
{
};


struct DummyVecField
{
    static constexpr std::size_t dimension = 1;
    using field_type                       = DummyField;
    DummyVecField(std::string name, HybridQuantity::Vector v) {}
    bool isUsable() const { return false; }
    bool isSettable() const { return true; }
};


struct DummyParticleInitializer
{
};


struct AnIonPopulation : public ::testing::Test
{
    IonPopulation<ParticleArray<1>, DummyVecField, DummyParticleInitializer> protons{"protons", 1.,
                                                                                     nullptr};
    virtual ~AnIonPopulation();
};

AnIonPopulation::~AnIonPopulation() {}



TEST_F(AnIonPopulation, hasAMass)
{
    EXPECT_DOUBLE_EQ(1., protons.mass());
}




TEST_F(AnIonPopulation, hasAName)
{
    EXPECT_EQ("protons", protons.name());
}




TEST_F(AnIonPopulation, isNonUsableUponConstruction)
{
    EXPECT_FALSE(protons.isUsable());
}




TEST_F(AnIonPopulation, isSettableIfNonUsable)
{
    if (!protons.isUsable())
    {
        EXPECT_TRUE(protons.isSettable());
    }
}




TEST_F(AnIonPopulation, throwsIfOneWantsToAccessParticleBuffersWhileNotUsable)
{
    EXPECT_ANY_THROW(protons.domainParticles());
    EXPECT_ANY_THROW(protons.ghostParticles());
    EXPECT_ANY_THROW(protons.coarseToFineParticles());
}




TEST_F(AnIonPopulation, isResourceUserAndHasGetParticleArrayNamesOK)
{
    auto bufferNames = protons.getParticleArrayNames();
    EXPECT_EQ(1, bufferNames.size());
    EXPECT_EQ(protons.name(), bufferNames[0].name);
}



TEST_F(AnIonPopulation, isResourceUserAndHasFieldNamesAndQuantitiesOK)
{
    auto fieldProperties = protons.getFieldNamesAndQuantities();
    EXPECT_EQ(protons.name() + std::string{"_rho"}, fieldProperties[0].name);
    EXPECT_EQ(HybridQuantity::Scalar::rho, fieldProperties[0].qty);
}



TEST_F(AnIonPopulation, hasAVecFieldSubResource)
{
    DummyVecField const& vf = std::get<0>(protons.getCompileTimeResourcesUserList());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
