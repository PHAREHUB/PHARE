
#include "amr/types/amr_types.h"
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/maxwellian_particle_initializer.h"
#include "data/vecfield/vecfield.h"
#include "data_provider.h"
#include "messengers/hybrid_messenger.h"
#include "messengers/messenger.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "resources_manager/resources_manager.h"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;
using namespace PHARE::solver;
using namespace PHARE::amr;


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using VecField1D                         = VecField<NdArrayVector1D<>, HybridQuantity>;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using MaxwellianParticleInitializer1D
    = MaxwellianParticleInitializer<ParticleArray<dim>, GridYee1D>;
using IonsPop1D         = IonPopulation<ParticleArray<dim>, VecField1D, GridYee1D>;
using Ions1D            = Ions<IonsPop1D, GridYee1D>;
using Electromag1D      = Electromag<VecField1D>;
using HybridModelT      = HybridModel<GridYee1D, Electromag1D, Ions1D, SAMRAI_Types>;
using MHDModelT         = MHDModel<GridYee1D, VecField1D, SAMRAI_Types>;
using ResourcesManagerT = ResourcesManager<GridYee1D>;




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

PHARE::initializer::PHAREDict createIonsDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["ions"]["name"]           = std::string{"ions"};
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

    return dict;
}




TEST(AHybridModel, fillsHybridMessengerInfo)
{
    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::unique_ptr<HybridModelT> hybridModel{
        std::make_unique<HybridModelT>(createIonsDict(), resourcesManagerHybrid)};




    std::unique_ptr<IMessengerInfo> modelInfoPtr = std::make_unique<HybridMessengerInfo>();

    hybridModel->fillMessengerInfo(modelInfoPtr);

    auto& modelInfo = dynamic_cast<HybridMessengerInfo const&>(*modelInfoPtr);


    EXPECT_EQ("EM_B", modelInfo.modelMagnetic.vecName);
    EXPECT_EQ("EM_B_x", modelInfo.modelMagnetic.xName);
    EXPECT_EQ("EM_B_y", modelInfo.modelMagnetic.yName);
    EXPECT_EQ("EM_B_z", modelInfo.modelMagnetic.zName);

    EXPECT_EQ("EM_E", modelInfo.modelElectric.vecName);
    EXPECT_EQ("EM_E_x", modelInfo.modelElectric.xName);
    EXPECT_EQ("EM_E_y", modelInfo.modelElectric.yName);
    EXPECT_EQ("EM_E_z", modelInfo.modelElectric.zName);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();

    // Finalize
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
