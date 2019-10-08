

#include <string>

#include "data/grid/gridlayoutdefs.h"




#include "utilities/index/index.h"


#include "data_provider.h"
#include "python_data_provider.h"
#include "restart_data_provider.h"




#include "data/electromag/electromag.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutimplyee.h"
#include "data/particles/particle_array.h"

#include <memory>

using namespace PHARE::initializer;

using GridLayoutT    = PHARE::core::GridLayout<PHARE::core::GridLayoutImplYee<1, 1>>;
using ParticleArrayT = PHARE::core::ParticleArray<1>;

#include "gmock/gmock.h"
#include "gtest/gtest.h"




TEST(APythonDataProvider, isADataProvider)
{
    PHAREDictHandler::INSTANCE().init<1>();
    char const* argv                   = "job_super.py";
    std::unique_ptr<DataProvider> pydp = std::make_unique<PythonDataProvider>(2, argv);
    PHAREDictHandler::INSTANCE().stop<1>();
}


TEST(ARestartTimeDataProvider, isADataProvider)
{
    std::unique_ptr<DataProvider> rdp = std::make_unique<RestartDataProvider>();
}



TEST(APythonDataProvider, providesAValidTree)
{
    PHAREDictHandler::INSTANCE().init<1>();

    char const* name                         = "init.py";
    std::unique_ptr<PythonDataProvider> pydp = std::make_unique<PythonDataProvider>(2, name);
    pydp->read();
    auto& input = PHAREDictHandler::INSTANCE().dict<1>();

    auto simulationName = input["simulation"]["name"].to<std::string>();
    auto nbrPopulations = input["simulation"]["ions"]["nbr_populations"].to<int>();

    auto& pop0                       = input["simulation"]["ions"]["pop0"];
    auto pop0Name                    = pop0["name"].to<std::string>();
    auto pop0Mass                    = pop0["mass"].to<double>();
    auto& pop0ParticleInitializer    = pop0["particle_initializer"];
    auto pop0ParticleInitializerName = pop0ParticleInitializer["name"].to<std::string>();
    auto pop0density                 = pop0ParticleInitializer["density"].to<ScalarFunction<1>>();
    auto bulk0x             = pop0ParticleInitializer["bulk_velocity_x"].to<ScalarFunction<1>>();
    auto bulk0y             = pop0ParticleInitializer["bulk_velocity_y"].to<ScalarFunction<1>>();
    auto bulk0z             = pop0ParticleInitializer["bulk_velocity_z"].to<ScalarFunction<1>>();
    auto vth0x              = pop0ParticleInitializer["thermal_velocity_x"].to<ScalarFunction<1>>();
    auto vth0y              = pop0ParticleInitializer["thermal_velocity_y"].to<ScalarFunction<1>>();
    auto vth0z              = pop0ParticleInitializer["thermal_velocity_z"].to<ScalarFunction<1>>();
    auto pop0NbrPartPerCell = pop0ParticleInitializer["nbr_part_per_cell"].to<int>();
    auto pop0Charge         = pop0ParticleInitializer["charge"].to<double>();
    auto pop0Basis          = pop0ParticleInitializer["basis"].to<std::string>();


    EXPECT_EQ("simulation_test", simulationName);
    EXPECT_EQ(2, nbrPopulations);
    EXPECT_EQ("protons", pop0Name);
    EXPECT_DOUBLE_EQ(1., pop0Mass);
    EXPECT_EQ("maxwellian", pop0ParticleInitializerName);
    EXPECT_DOUBLE_EQ(4., pop0density(2.));
    EXPECT_DOUBLE_EQ(4., bulk0x(2.));
    EXPECT_DOUBLE_EQ(6., bulk0y(2.));
    EXPECT_DOUBLE_EQ(8., bulk0z(2.));
    EXPECT_DOUBLE_EQ(10., vth0x(2.));
    EXPECT_DOUBLE_EQ(12., vth0y(2.));
    EXPECT_DOUBLE_EQ(14., vth0z(2.));
    EXPECT_EQ(100, pop0NbrPartPerCell);
    EXPECT_DOUBLE_EQ(1., pop0Charge);
    EXPECT_EQ("cartesian", pop0Basis);

    PHAREDictHandler::INSTANCE().stop<1>();
}



int main(int argc, char** argv)

{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
