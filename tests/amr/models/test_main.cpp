
#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/ions/ion_initializer.h"
#include "data/ions/ion_population/ion_population.h"
#include "data/ions/ions.h"
#include "data/ions/particle_initializers/fluid_particle_initializer.h"
#include "data/vecfield/vecfield.h"
#include "evolution/messengers/hybrid_messenger.h"
#include "evolution/messengers/messenger.h"
#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "tools/resources_manager.h"


#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE;


static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;
using VecField1D                         = VecField<NdArrayVector1D<>, HybridQuantity>;
using IonsPop1D                          = IonPopulation<ParticleArray<dim>, VecField1D>;
using GridImplYee1D                      = GridLayoutImplYee<dim, interpOrder>;
using GridYee1D                          = GridLayout<GridImplYee1D>;
using Ions1D                             = Ions<IonsPop1D, GridYee1D>;
using Electromag1D                       = Electromag<VecField1D>;
using IonsInit1D                         = IonsInitializer<ParticleArray<dim>, GridYee1D>;
using FluidParticleInitializer1D         = FluidParticleInitializer<ParticleArray<dim>, GridYee1D>;
using HybridModelT                       = HybridModel<GridYee1D, Electromag1D, Ions1D, IonsInit1D>;
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




auto getIonsInit()
{
    IonsInit1D ionsInit;
    ionsInit.name   = "Ions";
    ionsInit.masses = {{0.1, 0.3}};

    ionsInit.names.emplace_back("specie1");
    ionsInit.names.emplace_back("specie2");

    ionsInit.nbrPopulations = 2;

    ionsInit.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
        std::make_unique<ScalarFunction<dim>>(density),
        std::make_unique<VectorFunction<dim>>(bulkVelocity),
        std::make_unique<VectorFunction<dim>>(thermalVelocity), -1., 10));

    ionsInit.particleInitializers.push_back(std::make_unique<FluidParticleInitializer1D>(
        std::make_unique<ScalarFunction<dim>>(density),
        std::make_unique<VectorFunction<dim>>(bulkVelocity),
        std::make_unique<VectorFunction<dim>>(thermalVelocity), -1., 10));

    return ionsInit;
}




TEST(AHybridModel, fillsHybridMessengerInfo)
{
    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::unique_ptr<HybridModelT> hybridModel{
        std::make_unique<HybridModelT>(getIonsInit(), resourcesManagerHybrid)};




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

    EXPECT_EQ("Ions_rho", modelInfo.modelIonDensity);
    EXPECT_EQ("Ions_bulkVel", modelInfo.modelIonBulk.vecName);
    EXPECT_EQ("Ions_bulkVel_x", modelInfo.modelIonBulk.xName);
    EXPECT_EQ("Ions_bulkVel_y", modelInfo.modelIonBulk.yName);
    EXPECT_EQ("Ions_bulkVel_z", modelInfo.modelIonBulk.zName);



    EXPECT_NE(std::end(modelInfo.initIonDensity),
              std::find(std::begin(modelInfo.initIonDensity), std::end(modelInfo.initIonDensity),
                        "Ions_rho"));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.vecName == std::string{"Ions_bulkVel"}; }));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.xName == std::string{"Ions_bulkVel_x"}; }));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.yName == std::string{"Ions_bulkVel_y"}; }));

    EXPECT_NE(
        std::end(modelInfo.initIonBulk),
        std::find_if(std::begin(modelInfo.initIonBulk), std::end(modelInfo.initIonBulk),
                     [](auto const& desc) { return desc.zName == std::string{"Ions_bulkVel_z"}; }));
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
