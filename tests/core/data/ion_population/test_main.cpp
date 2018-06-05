
#include <type_traits>


#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "data/ion_population/ion_population.h"
#include "hybrid/hybrid_quantities.h"

using namespace PHARE;




struct DummyParticleArray
{
};


struct DummyField
{
};


struct DummyVecField
{
    static constexpr std::size_t dimension = 1;
    using field_type                       = DummyField;
    DummyVecField(std::string name, HybridQuantity::Vector v) {}
};



struct AnIonPopulation : public ::testing::Test
{
    IonPopulation<DummyParticleArray, DummyVecField> protons{"protons", 1.};
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
    EXPECT_EQ(3, bufferNames.size());
    EXPECT_EQ(protons.name() + std::string{"_domain"}, bufferNames[0].name);
    EXPECT_EQ(protons.name() + std::string{"_ghost"}, bufferNames[1].name);
    EXPECT_EQ(protons.name() + std::string{"_coarseToFine"}, bufferNames[2].name);
}



TEST_F(AnIonPopulation, isResourceUserAndHasFieldNamesAndQuantitiesOK)
{
    auto fieldProperties = protons.getFieldNamesAndQuantities();
    EXPECT_EQ(protons.name() + std::string{"_rho"}, fieldProperties[0].name);
    EXPECT_EQ(HybridQuantity::Scalar::rho, fieldProperties[0].qty);
}



TEST_F(AnIonPopulation, hasAVecFieldSubResource)
{
    DummyVecField const& vf = std::get<0>(protons.getSubResourcesObject());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
