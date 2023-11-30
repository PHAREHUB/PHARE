#include <type_traits>




#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/particles/particle_array.hpp"
#include "initializer/data_provider.hpp"
#include "core/hybrid/hybrid_quantities.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;
using namespace PHARE::initializer;



struct DummyField
{
};


struct DummyVecField
{
    static constexpr std::size_t dimension = 1;
    using grid_type                        = DummyField;
    using field_type                       = DummyField;
    DummyVecField(std::string name, [[maybe_unused]] HybridQuantity::Vector v) { (void)name; }
    bool isUsable() const { return false; }
    bool isSettable() const { return true; }
};


struct DummyParticleInitializer
{
};


struct DummyLayout
{
};

PHAREDict getDict()
{
    PHAREDict dict;
    dict["name"]                         = std::string{"protons"};
    dict["mass"]                         = 1.;
    dict["particle_initializer"]["name"] = std::string{"DummyParticleInitializer"};
    return dict;
}

struct AnIonPopulation : public ::testing::Test
{
    IonPopulation<ParticleArray<1>, DummyVecField, DummyLayout> protons{getDict()};
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
    EXPECT_ANY_THROW(auto& p = protons.domainParticles());
    EXPECT_ANY_THROW(auto& p = protons.patchGhostParticles());
    EXPECT_ANY_THROW(auto& p = protons.levelGhostParticles());
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
    [[maybe_unused]] DummyVecField const& vf
        = std::get<0>(protons.getCompileTimeResourcesUserList());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
