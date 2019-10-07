

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
    auto& input = PHAREDictHandler::INSTANCE().dict<1>(); // dict<1>();

    EXPECT_EQ("simulation_test", input["simulation"]["name"].to<std::string>());
    EXPECT_EQ(2, input["simulation"]["ions"]["nbr_populations"].to<int>());

    auto& pop0 = input["simulation"]["ions"]["pop0"];
    EXPECT_EQ("protons", pop0["name"].to<std::string>());
    EXPECT_EQ(1., pop0["mass"].to<double>());
    EXPECT_EQ(1., pop0["charge"].to<double>());
    EXPECT_EQ("maxwellian", pop0["particle_initializer"]["name"].to<std::string>());

    auto density0
        = pop0["particle_initializer"]["density"].to<PHARE::initializer::ScalarFunction<1>>();
    auto d = density0(2.);
    EXPECT_EQ(4., d);
    PHAREDictHandler::INSTANCE().stop<1>();
}



int main(int argc, char** argv)

{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
