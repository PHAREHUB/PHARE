#include "amr/types/amr_types.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ions/ion_population/ion_population.hpp"
#include "core/data/ions/ions.hpp"
#include "core/data/electrons/electrons.hpp"
#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/tensorfield/tensorfield.hpp"
#include "initializer/data_provider.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/resources_manager/resources_manager.hpp"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using namespace PHARE::core;
using namespace PHARE::solver;
using namespace PHARE::amr;

#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

static constexpr std::size_t dim         = 1;
static constexpr std::size_t interpOrder = 1;

using Field_t          = Field<1, HybridQuantity::Scalar>;
using Grid1D           = Grid<NdArrayVector<1>, HybridQuantity::Scalar>;
using VecField1D       = VecField<Field_t, HybridQuantity>;
using SymTensorField1D = SymTensorField<Field_t, HybridQuantity>;
using GridImplYee1D    = GridLayoutImplYee<dim, interpOrder>;
using ParticleArray1D  = ParticleArray<dim>;
using GridYee1D        = GridLayout<GridImplYee1D>;

using MaxwellianParticleInitializer1D = MaxwellianParticleInitializer<ParticleArray1D, GridYee1D>;

using IonsPop1D    = IonPopulation<ParticleArray1D, VecField1D, SymTensorField1D>;
using Ions1D       = Ions<IonsPop1D, GridYee1D>;
using Electromag1D = Electromag<VecField1D>;
using Electrons1D  = Electrons<Ions1D>;
using HybridModelT
    = HybridModel<GridYee1D, Electromag1D, Ions1D, Electrons1D, SAMRAI_Types, Grid1D>;
using MHDModelT         = MHDModel<GridYee1D, VecField1D, SAMRAI_Types, Grid1D>;
using ResourcesManagerT = ResourcesManager<GridYee1D, Grid1D>;


using InitFunctionT = PHARE::initializer::InitFunction<1>;

PHARE::initializer::PHAREDict createDict()
{
    PHARE::initializer::PHAREDict dict;
    dict["ions"]["nbrPopulations"] = std::size_t{2};
    dict["ions"]["pop0"]["name"]   = std::string{"protons"};
    dict["ions"]["pop0"]["mass"]   = 1.;
    dict["ions"]["pop0"]["particle_initializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop0"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(vx);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(vy);

    dict["ions"]["pop0"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(vz);


    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(vthx);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(vthy);

    dict["ions"]["pop0"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(vthz);


    dict["ions"]["pop0"]["particle_initializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop0"]["particle_initializer"]["charge"]         = -1.;
    dict["ions"]["pop0"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

    dict["ions"]["pop1"]["name"] = std::string{"alpha"};
    dict["ions"]["pop1"]["mass"] = 1.;
    dict["ions"]["pop1"]["particle_initializer"]["name"]
        = std::string{"MaxwellianParticleInitializer"};
    dict["ions"]["pop1"]["particle_initializer"]["density"] = static_cast<InitFunctionT>(density);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_x"]
        = static_cast<InitFunctionT>(vx);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_y"]
        = static_cast<InitFunctionT>(vy);

    dict["ions"]["pop1"]["particle_initializer"]["bulk_velocity_z"]
        = static_cast<InitFunctionT>(vz);


    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_x"]
        = static_cast<InitFunctionT>(vthx);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_y"]
        = static_cast<InitFunctionT>(vthy);

    dict["ions"]["pop1"]["particle_initializer"]["thermal_velocity_z"]
        = static_cast<InitFunctionT>(vthz);


    dict["ions"]["pop1"]["particle_initializer"]["nbrPartPerCell"] = int{100};
    dict["ions"]["pop1"]["particle_initializer"]["charge"]         = -1.;
    dict["ions"]["pop1"]["particle_initializer"]["basis"]          = std::string{"Cartesian"};

    dict["electromag"]["name"]             = std::string{"EM"};
    dict["electromag"]["electric"]["name"] = std::string{"E"};
    dict["electromag"]["magnetic"]["name"] = std::string{"B"};

    dict["electromag"]["magnetic"]["initializer"]["x_component"] = static_cast<InitFunctionT>(bx);
    dict["electromag"]["magnetic"]["initializer"]["y_component"] = static_cast<InitFunctionT>(by);
    dict["electromag"]["magnetic"]["initializer"]["z_component"] = static_cast<InitFunctionT>(bz);

    dict["electrons"]["pressure_closure"]["name"] = std::string{"isothermal"};
    dict["electrons"]["pressure_closure"]["Te"]   = 0.12;

    return dict;
}




TEST(AHybridModel, fillsHybridMessengerInfo)
{
    std::shared_ptr<ResourcesManagerT> resourcesManagerHybrid{
        std::make_shared<ResourcesManagerT>()};

    std::unique_ptr<HybridModelT> hybridModel{
        std::make_unique<HybridModelT>(createDict(), resourcesManagerHybrid)};




    std::unique_ptr<IMessengerInfo> modelInfoPtr = std::make_unique<HybridMessengerInfo>();

    hybridModel->fillMessengerInfo(modelInfoPtr);

    auto& modelInfo = dynamic_cast<HybridMessengerInfo const&>(*modelInfoPtr);


    EXPECT_EQ("EM_B", modelInfo.modelMagnetic);
    // EXPECT_EQ("EM_B_x", modelInfo.modelMagnetic.xName);
    // EXPECT_EQ("EM_B_y", modelInfo.modelMagnetic.yName);
    // EXPECT_EQ("EM_B_z", modelInfo.modelMagnetic.zName);

    EXPECT_EQ("EM_E", modelInfo.modelElectric);
    // EXPECT_EQ("EM_E_x", modelInfo.modelElectric.xName);
    // EXPECT_EQ("EM_E_y", modelInfo.modelElectric.yName);
    // EXPECT_EQ("EM_E_z", modelInfo.modelElectric.zName);
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
