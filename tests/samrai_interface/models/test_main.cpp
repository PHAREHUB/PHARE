
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/fluid_particle_initializer.h"
#include "data/vecfield/vecfield.h"
#include "data_provider.h"
#include "evolution/messengers/hybrid_messenger.h"
#include "evolution/messengers/messenger.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "tools/resources_manager.h"


#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;
using namespace PHARE::amr_interface;


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using VecField1D                         = VecField<NdArrayVector1D<>, HybridQuantity>;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using FluidParticleInitializer1D         = FluidParticleInitializer<ParticleArray<dim>, GridYee1D>;
using IonsPop1D                          = IonPopulation<ParticleArray<dim>, VecField1D, GridYee1D>;
using Ions1D                             = Ions<IonsPop1D, GridYee1D>;
using Electromag1D                       = Electromag<VecField1D>;
using HybridModelT                       = HybridModel<GridYee1D, Electromag1D, Ions1D>;
using MHDModelT                          = MHDModel<GridYee1D, VecField1D>;
using ResourcesManagerT                  = ResourcesManager<GridYee1D>;




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




using ScalarFunction = PHARE::initializer::ScalarFunction<1>;
using VectorFunction = PHARE::initializer::VectorFunction<1>;

PHARE::initializer::PHAREDict<1> createIonsDict()
{
    PHARE::initializer::PHAREDict<1> dict;
    dict["ions"]["name"]                                = std::string{"ions"};
    dict["ions"]["nbrPopulations"]                      = std::size_t{2};
    dict["ions"]["pop0"]["name"]                        = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]                        = 1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["name"] = std::string{"FluidParticleInitializer"};
    dict["ions"]["pop0"]["ParticleInitializer"]["density"] = static_cast<ScalarFunction>(density);

    dict["ions"]["pop0"]["ParticleInitializer"]["bulkVelocity"]
        = static_cast<VectorFunction>(bulkVelocity);

    dict["ions"]["pop0"]["ParticleInitializer"]["thermalVelocity"]
        = static_cast<VectorFunction>(thermalVelocity);

    dict["ions"]["pop0"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
    dict["ions"]["pop0"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop0"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};



    dict["ions"]["pop1"]["name"]                        = std::string{"protons"};
    dict["ions"]["pop1"]["mass"]                        = 1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["name"] = std::string{"FluidParticleInitializer"};
    dict["ions"]["pop1"]["ParticleInitializer"]["density"] = static_cast<ScalarFunction>(density);

    dict["ions"]["pop1"]["ParticleInitializer"]["bulkVelocity"]
        = static_cast<VectorFunction>(bulkVelocity);

    dict["ions"]["pop1"]["ParticleInitializer"]["thermalVelocity"]
        = static_cast<VectorFunction>(thermalVelocity);

    dict["ions"]["pop1"]["ParticleInitializer"]["nbrPartPerCell"] = std::size_t{100};
    dict["ions"]["pop1"]["ParticleInitializer"]["charge"]         = -1.;
    dict["ions"]["pop1"]["ParticleInitializer"]["basis"]          = std::string{"Cartesian"};

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

    EXPECT_EQ("ions_rho", modelInfo.modelIonDensity);
    EXPECT_EQ("ions_bulkVel", modelInfo.modelIonBulk.vecName);
    EXPECT_EQ("ions_bulkVel_x", modelInfo.modelIonBulk.xName);
    EXPECT_EQ("ions_bulkVel_y", modelInfo.modelIonBulk.yName);
    EXPECT_EQ("ions_bulkVel_z", modelInfo.modelIonBulk.zName);



    EXPECT_NE(std::end(modelInfo.initIonDensity),
              std::find(std::begin(modelInfo.initIonDensity), std::end(modelInfo.initIonDensity),
                        "ions_rho"));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.vecName == std::string{"ions_bulkVel"}; }));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.xName == std::string{"ions_bulkVel_x"}; }));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.yName == std::string{"ions_bulkVel_y"}; }));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.zName == std::string{"ions_bulkVel_z"}; }));
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
