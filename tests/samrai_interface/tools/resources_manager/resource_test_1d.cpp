
#include "resource_test_1d.h"
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data_provider.h"
#include "models/hybrid_state.h"


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;

using VecField1D                      = VecField<NdArrayVector1D<>, HybridQuantity>;
using IonPopulation1D                 = IonPopulation<ParticleArray<1>, VecField1D, GridYee1D>;
using Ions1D                          = Ions<IonPopulation1D, GridYee1D>;
using Electromag1D                    = Electromag<VecField1D>;
using MaxwellianParticleInitializer1D = MaxwellianParticleInitializer<ParticleArray<1>, GridYee1D>;
using HybridState1D                   = HybridState<Electromag1D, Ions1D>;


PHARE::initializer::PHAREDict<1> getDict()
{
    PHARE::initializer::PHAREDict<1> dict;
    dict["name"]                        = std::string{"protons"};
    dict["mass"]                        = 1.;
    dict["ParticleInitializer"]["name"] = std::string{"DummyParticleInitializer"};
    return dict;
}



struct IonPopulation1D_P
{
    std::string ionName = "ions";
    IonPopulation1D user{ionName, getDict()};
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
    (void)x;
    return std::array<double, 3>{{1.0, 0.0, 0.0}};
}


std::array<double, 3> thermalVelocity(double x)
{
    (void)x;
    return std::array<double, 3>{{0.5, 0.0, 0.0}};
}




struct Ions1D_P
{
    using ScalarFunction = PHARE::initializer::ScalarFunction<1>;
    using VectorFunction = PHARE::initializer::VectorFunction<1>;

    PHARE::initializer::PHAREDict<1> createIonsDict()
    {
        PHARE::initializer::PHAREDict<1> dict;
        dict["name"]                                = std::string{"ions"};
        dict["nbrPopulations"]                      = std::size_t{2};
        dict["pop0"]["name"]                        = std::string{"protons"};
        dict["pop0"]["mass"]                        = 1.;
        dict["pop0"]["ParticleInitializer"]["name"] = std::string{"MaxwellianParticleInitializer"};
        dict["pop0"]["ParticleInitializer"]["density"] = static_cast<ScalarFunction>(density);

        dict["pop0"]["ParticleInitializer"]["bulkVelocity"]
            = static_cast<VectorFunction>(bulkVelocity);

        dict["pop0"]["ParticleInitializer"]["thermalVelocity"]
            = static_cast<VectorFunction>(thermalVelocity);

        dict["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
        dict["pop0"]["ParticleInitializer"]["charge"]         = -1.;
        dict["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};



        dict["pop1"]["name"]                        = std::string{"protons"};
        dict["pop1"]["mass"]                        = 1.;
        dict["pop1"]["ParticleInitializer"]["name"] = std::string{"MaxwellianParticleInitializer"};
        dict["pop1"]["ParticleInitializer"]["density"] = static_cast<ScalarFunction>(density);

        dict["pop1"]["ParticleInitializer"]["bulkVelocity"]
            = static_cast<VectorFunction>(bulkVelocity);

        dict["pop1"]["ParticleInitializer"]["thermalVelocity"]
            = static_cast<VectorFunction>(thermalVelocity);

        dict["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
        dict["pop1"]["ParticleInitializer"]["charge"]         = -1.;
        dict["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

        return dict;
    }

    Ions1D user{createIonsDict()};
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
    using ScalarFunction = PHARE::initializer::ScalarFunction<1>;
    using VectorFunction = PHARE::initializer::VectorFunction<1>;

    PHARE::initializer::PHAREDict<1> createIonsDict()
    {
        PHARE::initializer::PHAREDict<1> dict;
        dict["ions"]["name"]           = std::string{"ions"};
        dict["ions"]["nbrPopulations"] = std::size_t{2};
        dict["ions"]["pop0"]["name"]   = std::string{"protons"};
        dict["ions"]["pop0"]["mass"]   = 1.;
        dict["ions"]["pop0"]["ParticleInitializer"]["name"]
            = std::string{"MaxwellianParticleInitializer"};
        dict["ions"]["pop0"]["ParticleInitializer"]["density"]
            = static_cast<ScalarFunction>(density);

        dict["ions"]["pop0"]["ParticleInitializer"]["bulkVelocity"]
            = static_cast<VectorFunction>(bulkVelocity);

        dict["ions"]["pop0"]["ParticleInitializer"]["thermalVelocity"]
            = static_cast<VectorFunction>(thermalVelocity);

        dict["ions"]["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
        dict["ions"]["pop0"]["ParticleInitializer"]["charge"]         = -1.;
        dict["ions"]["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};



        dict["ions"]["pop1"]["name"] = std::string{"protons"};
        dict["ions"]["pop1"]["mass"] = 1.;
        dict["ions"]["pop1"]["ParticleInitializer"]["name"]
            = std::string{"MaxwellianParticleInitializer"};
        dict["ions"]["pop1"]["ParticleInitializer"]["density"]
            = static_cast<ScalarFunction>(density);

        dict["ions"]["pop1"]["ParticleInitializer"]["bulkVelocity"]
            = static_cast<VectorFunction>(bulkVelocity);

        dict["ions"]["pop1"]["ParticleInitializer"]["thermalVelocity"]
            = static_cast<VectorFunction>(thermalVelocity);

        dict["ions"]["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
        dict["ions"]["pop1"]["ParticleInitializer"]["charge"]         = -1.;
        dict["ions"]["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

        return dict;
    }


    HybridState1D user{createIonsDict()};
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




TEST(usingResourcesManager, toGetTimeOfAResourcesUser)
{
    std::unique_ptr<BasicHierarchy> hierarchy;
    ResourcesManager<GridLayout<GridLayoutImplYee<1, 1>>> resourcesManager;
    IonPopulation1D_P pop;
    auto s    = inputBase + std::string("/input/input_db_1d");
    hierarchy = std::make_unique<BasicHierarchy>(inputBase + std::string("/input/input_db_1d"));
    hierarchy->init();
    resourcesManager.registerResources(pop.user);
    auto& patchHierarchy = hierarchy->hierarchy;

    double const initDataTime{3.14};

    for (int iLevel = 0; iLevel < patchHierarchy->getNumberOfLevels(); ++iLevel)
    {
        auto patchLevel = patchHierarchy->getPatchLevel(iLevel);
        for (auto& patch : *patchLevel)
        {
            resourcesManager.allocate(pop.user, *patch, initDataTime);
            auto times = resourcesManager.getTime(pop.user, *patch);

            EXPECT_TRUE(std::equal(std::begin(times) + 1, std::end(times), std::begin(times)));
            EXPECT_DOUBLE_EQ(initDataTime, times[0]);
        }
    }
}




REGISTER_TYPED_TEST_CASE_P(aResourceUserCollection, hasPointersValidOnlyWithGuard);


typedef ::testing::Types<IonPop1DOnly, VecField1DOnly, Ions1DOnly, Electromag1DOnly,
                         HybridState1DOnly>
    MyTypes;
INSTANTIATE_TYPED_TEST_CASE_P(testResourcesManager, aResourceUserCollection, MyTypes);
