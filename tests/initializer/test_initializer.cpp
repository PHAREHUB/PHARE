
#include <memory>
#include <string>
#include <stdlib.h>


#include "initializer/data_provider.hpp"
#include "initializer/python_data_provider.hpp"
#include "initializer/restart_data_provider.hpp"


#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/electromag/electromag.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/data/particles/particle_array.hpp"



using namespace PHARE::initializer;

using GridLayoutT    = PHARE::core::GridLayout<PHARE::core::GridLayoutImplYee<1, 1>>;
using ParticleArrayT = PHARE::core::ParticleArray<1>;

#include "gmock/gmock.h"
#include "gtest/gtest.h"




TEST(APythonDataProvider, isADataProvider)
{
    PHAREDictHandler::INSTANCE().init();
    char const* argv                   = "job";
    std::unique_ptr<DataProvider> pydp = std::make_unique<PythonDataProvider>(argv);
    PHAREDictHandler::INSTANCE().stop();
}


TEST(ARestartTimeDataProvider, isADataProvider)
{
    std::unique_ptr<DataProvider> rdp = std::make_unique<RestartDataProvider>();
}



TEST(APythonDataProvider, providesAValidTree)
{
    PHAREDictHandler::INSTANCE().init();

    char const* name = "job";
    PythonDataProvider pydp{name};
    pydp.read();
    auto& input = PHAREDictHandler::INSTANCE().dict();

    auto simulationName = input["simulation"]["name"].to<std::string>();
    auto dim            = input["simulation"]["dimension"].to<int>();
    auto interp_order   = input["simulation"]["interp_order"].to<int>();
    auto dt             = input["simulation"]["time_step"].to<double>();

    auto layout = input["simulation"]["grid"]["layout_type"].to<std::string>();
    auto nx     = input["simulation"]["grid"]["nbr_cells"]["x"].to<int>();
    auto dx     = input["simulation"]["grid"]["meshsize"]["x"].to<double>();
    auto origin = input["simulation"]["grid"]["origin"]["x"].to<double>();

    auto pusherName
        = input["simulation"]["algo"]["ion_updater"]["pusher"]["name"].to<std::string>();


    auto nbrPopulations              = input["simulation"]["ions"]["nbrPopulations"].to<int>();
    auto& pop0                       = input["simulation"]["ions"]["pop0"];
    auto pop0Name                    = pop0["name"].to<std::string>();
    auto pop0Mass                    = pop0["mass"].to<double>();
    auto& pop0ParticleInitializer    = pop0["particle_initializer"];
    auto pop0ParticleInitializerName = pop0ParticleInitializer["name"].to<std::string>();
    auto pop0density                 = pop0ParticleInitializer["density"].to<InitFunction<1>>();
    auto bulk0x             = pop0ParticleInitializer["bulk_velocity_x"].to<InitFunction<1>>();
    auto bulk0y             = pop0ParticleInitializer["bulk_velocity_y"].to<InitFunction<1>>();
    auto bulk0z             = pop0ParticleInitializer["bulk_velocity_z"].to<InitFunction<1>>();
    auto vth0x              = pop0ParticleInitializer["thermal_velocity_x"].to<InitFunction<1>>();
    auto vth0y              = pop0ParticleInitializer["thermal_velocity_y"].to<InitFunction<1>>();
    auto vth0z              = pop0ParticleInitializer["thermal_velocity_z"].to<InitFunction<1>>();
    auto pop0NbrPartPerCell = pop0ParticleInitializer["nbr_part_per_cell"].to<int>();
    auto pop0Charge         = pop0ParticleInitializer["charge"].to<double>();
    auto pop0Basis          = pop0ParticleInitializer["basis"].to<std::string>();


    EXPECT_EQ("simulation_test", simulationName);
    EXPECT_EQ(1, interp_order);

    EXPECT_EQ(1, dim);
    EXPECT_EQ(65, nx);
    EXPECT_DOUBLE_EQ(1. / 65., dx);
    EXPECT_DOUBLE_EQ(0.001, dt);
    EXPECT_DOUBLE_EQ(0., origin);
    EXPECT_EQ("yee", layout);

    EXPECT_EQ("modified_boris", pusherName);

    std::vector<double> input_2 = {2};

    EXPECT_EQ(2, nbrPopulations);
    EXPECT_EQ("protons", pop0Name);
    EXPECT_DOUBLE_EQ(1., pop0Mass);
    EXPECT_EQ("maxwellian", pop0ParticleInitializerName);
    EXPECT_DOUBLE_EQ(2., (*pop0density(input_2))[0]);
    EXPECT_DOUBLE_EQ(1., (*bulk0x(input_2))[0]);
    EXPECT_DOUBLE_EQ(1., (*bulk0y(input_2))[0]);
    EXPECT_DOUBLE_EQ(1., (*bulk0z(input_2))[0]);
    EXPECT_DOUBLE_EQ(1., (*vth0x(input_2))[0]);
    EXPECT_DOUBLE_EQ(1., (*vth0y(input_2))[0]);
    EXPECT_DOUBLE_EQ(1., (*vth0z(input_2))[0]);
    EXPECT_EQ(100, pop0NbrPartPerCell);
    EXPECT_DOUBLE_EQ(1., pop0Charge);
    EXPECT_EQ("cartesian", pop0Basis);

    PHAREDictHandler::INSTANCE().stop();
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
