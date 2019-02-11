
#include "resource_test_1d.h"
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/particle_initializers/fluid_particle_initializer.h"
#include "models/hybrid_state.h"


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;

using VecField1D                 = VecField<NdArrayVector1D<>, HybridQuantity>;
using IonPopulation1D            = IonPopulation<ParticleArray<1>, VecField1D, GridYee1D>;
using Ions1D                     = Ions<IonPopulation1D, GridYee1D>;
using Electromag1D               = Electromag<VecField1D>;
using FluidParticleInitializer1D = FluidParticleInitializer<ParticleArray<1>, GridYee1D>;
using IonInitializer1D           = IonsInitializer<ParticleArray<1>, GridYee1D>;
using HybridState1D              = HybridState<Electromag1D, Ions1D, IonInitializer1D>;


struct IonPopulation1D_P
{
    std::string name = "protons";
    double mass      = 1.;
    IonPopulation1D user{name, mass, nullptr};
};


struct VecField1D_P
{
    std::string name = "B";
    HybridQuantity::Vector qty{HybridQuantity::Vector::B};
    VecField1D user{name, qty};
};



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




struct Ions1D_P
{
    IonsInitializer<ParticleArray<1>, GridYee1D> createInitializer()
    {
        IonsInitializer<ParticleArray<1>, GridYee1D> initializer;

        initializer.masses.push_back(1.);
        initializer.names.push_back("protons");
        initializer.nbrPopulations = 1;
        initializer.name           = "TestIons";
        initializer.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
            &density, &bulkVelocity, &thermalVelocity, -1., 10));


        return initializer;
    }


    Ions1D user{createInitializer()};
};




struct Electromag1D_P
{
    std::string name = "ElectroTest";
    Electromag1D user;
    Electromag1D_P()
        : user{name}
    {
    }
};



struct HybridState1D_P
{
    IonsInitializer<ParticleArray<1>, GridYee1D> createInitializer()
    {
        IonsInitializer<ParticleArray<1>, GridYee1D> initializer;

        initializer.masses.push_back(1.);
        initializer.names.push_back("protons");
        initializer.nbrPopulations = 1;
        initializer.name           = "TestIons";
        initializer.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
            &density, &bulkVelocity, &thermalVelocity, -1., 10));


        return initializer;
    }


    HybridState1D user{createInitializer()};
};




using IonPop1DOnly          = std::tuple<IonPopulation1D_P>;
using VecField1DOnly        = std::tuple<VecField1D_P>;
using Ions1DOnly            = std::tuple<Ions1D_P>;
using VecField1DAndIonPop1D = std::tuple<VecField1D_P, IonPopulation1D_P>;
using Electromag1DOnly      = std::tuple<Electromag1D_P>;
using HybridState1DOnly     = std::tuple<HybridState1D_P>;

TYPED_TEST_CASE_P(aResourceUserCollection);




TYPED_TEST_P(aResourceUserCollection, hasPointersValidOnlyWithGuard)
{
    TypeParam resourceUserCollection;

    auto check = [this](auto& resourceUserPack) {
        auto& hierarchy    = this->hierarchy->hierarchy;
        auto& resourceUser = resourceUserPack.user;

        for (int iLevel = 0; iLevel < hierarchy->getNumberOfLevels(); ++iLevel)
        {
            auto patchLevel = hierarchy->getPatchLevel(iLevel);
            for (auto const& patch : *patchLevel)
            {
                auto dataOnPatch = this->resourcesManager.setOnPatch(*patch, resourceUser);
                EXPECT_TRUE(resourceUser.isUsable());
                EXPECT_FALSE(resourceUser.isSettable());
            }
            EXPECT_FALSE(resourceUser.isUsable());
            EXPECT_TRUE(resourceUser.isSettable());
        }
    };

    std::apply(check, resourceUserCollection);
}



REGISTER_TYPED_TEST_CASE_P(aResourceUserCollection, hasPointersValidOnlyWithGuard);


typedef ::testing::Types<IonPop1DOnly, VecField1DOnly, Ions1DOnly, Electromag1DOnly,
                         HybridState1DOnly>
    MyTypes;
INSTANTIATE_TYPED_TEST_CASE_P(testResourcesManager, aResourceUserCollection, MyTypes);
