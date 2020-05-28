
#include "resource_test_1d.h"
#include "core/data/electromag/electromag.h"
#include "core/data/electrons/electrons.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "initializer/data_provider.h"
#include "core/models/hybrid_state.h"


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;

using VecField1D                      = VecField<NdArrayVector<1>, HybridQuantity>;
using IonPopulation1D                 = IonPopulation<ParticleArray<1>, VecField1D, GridYee1D>;
using Ions1D                          = Ions<IonPopulation1D, GridYee1D>;
using Electromag1D                    = Electromag<VecField1D>;
using Electrons1D                     = Electrons<Ions1D>;
using MaxwellianParticleInitializer1D = MaxwellianParticleInitializer<ParticleArray<1>, GridYee1D>;
using HybridState1D                   = HybridState<Electromag1D, Ions1D, Electrons1D>;



double density(double x)
{
    return x * x + 2.;
}

double vx(double x)
{
    (void)x;
    return 1.;
}


double vy(double x)
{
    (void)x;
    return 1.;
}


double vz(double x)
{
    (void)x;
    return 1.;
}


double vthx(double x)
{
    (void)x;
    return 1.;
}


double vthy(double x)
{
    (void)x;
    return 1.;
}



double vthz(double x)
{
    (void)x;
    return 1.;
}



double bx(double x)
{
    return x * x + 2.;
}

double by(double x)
{
    return x * x + 2.;
}

double bz(double x)
{
    return x * x + 2.;
}

double ex(double x)
{
    return x * x + 2.;
}

double ey(double x)
{
    return x * x + 2.;
}

double ez(double x)
{
    return x * x + 2.;
}




using ScalarFunctionT = PHARE::initializer::ScalarFunction<1>;

PHARE::initializer::PHAREDict createInitDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"] = int{2};
    dict["ions"]["pop0"]["name"]   = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]   = 1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop0"]["ParticleInitializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop0"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

    dict["ions"]["pop1"]["name"] = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"] = 1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop1"]["ParticleInitializer"]["density"] = static_cast<ScalarFunctionT>(density);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_x"]
        = static_cast<ScalarFunctionT>(vx);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_y"]
        = static_cast<ScalarFunctionT>(vy);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulk_velocity_z"]
        = static_cast<ScalarFunctionT>(vz);


    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_x"]
        = static_cast<ScalarFunctionT>(vthx);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_y"]
        = static_cast<ScalarFunctionT>(vthy);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermal_velocity_z"]
        = static_cast<ScalarFunctionT>(vthz);


    dict["ions"]["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop1"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["electric"]["initializer"]["x_component"] = static_cast<ScalarFunctionT>(ex);
    dict["electromag"]["electric"]["initializer"]["y_component"] = static_cast<ScalarFunctionT>(ey);
    dict["electromag"]["electric"]["initializer"]["z_component"] = static_cast<ScalarFunctionT>(ez);

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<ScalarFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<ScalarFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<ScalarFunctionT>(bz);

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = 0.12;


    return dict;
}


struct IonPopulation1D_P
{
    IonPopulation1D user{createInitDict()["ions"]["pop0"]};
};


struct VecField1D_P
{
    std::string name = "B";
    HybridQuantity::Vector qty{HybridQuantity::Vector::B};
    VecField1D user{name, qty};
};




struct Ions1D_P
{
    Ions1D user{createInitDict()["ions"]};
};




struct Electromag1D_P
{
    std::string name = "ElectroTest";
    Electromag1D user;
    Electromag1D_P()
        : user{createInitDict()["electromag"]}
    {
    }
};



struct HybridState1D_P
{
    HybridState1D user{createInitDict()};
};




using IonPop1DOnly          = std::tuple<IonPopulation1D_P>;
using VecField1DOnly        = std::tuple<VecField1D_P>;
using Ions1DOnly            = std::tuple<Ions1D_P>;
using VecField1DAndIonPop1D = std::tuple<VecField1D_P, IonPopulation1D_P>;
using Electromag1DOnly      = std::tuple<Electromag1D_P>;
using HybridState1DOnly     = std::tuple<HybridState1D_P>;

TYPED_TEST_SUITE_P(aResourceUserCollection);




TYPED_TEST_P(aResourceUserCollection, hasPointersValidOnlyWithGuard)
{
    TypeParam resourceUserCollection;

    auto check = [this](auto& resourceUserPack) {
        auto& hierarchy_   = this->hierarchy->hierarchy;
        auto& resourceUser = resourceUserPack.user;

        for (int iLevel = 0; iLevel < hierarchy_->getNumberOfLevels(); ++iLevel)
        {
            auto patchLevel = hierarchy_->getPatchLevel(iLevel);
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
            auto times = resourcesManager.getTimes(pop.user, *patch);

            EXPECT_TRUE(std::equal(std::begin(times) + 1, std::end(times), std::begin(times)));
            EXPECT_DOUBLE_EQ(initDataTime, times[0]);
        }
    }
}




REGISTER_TYPED_TEST_SUITE_P(aResourceUserCollection, hasPointersValidOnlyWithGuard);


typedef ::testing::Types<IonPop1DOnly, VecField1DOnly, Ions1DOnly, Electromag1DOnly,
                         HybridState1DOnly>
    MyTypes;
INSTANTIATE_TYPED_TEST_SUITE_P(testResourcesManager, aResourceUserCollection, MyTypes);
