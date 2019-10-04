

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

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>

using namespace PHARE::initializer;

using GridLayoutT    = PHARE::core::GridLayout<PHARE::core::GridLayoutImplYee<1, 1>>;
using ParticleArrayT = PHARE::core::ParticleArray<1>;



TEST(APythonDataProvider, isADataProvider)
{
    char const* argv                   = "job_super.py";
    std::unique_ptr<DataProvider> pydp = std::make_unique<PythonDataProvider>(2, argv);
}


TEST(ARestartTimeDataProvider, isADataProvider)
{
    std::unique_ptr<DataProvider> rdp = std::make_unique<RestartDataProvider>();
}



TEST(APythonDataProvider, providesAProperTree)
{
    char const* argv                         = "init.py";
    std::unique_ptr<PythonDataProvider> pydp = std::make_unique<PythonDataProvider>(2, argv);
    pydp->read();
    auto& input = dict<1>();

    EXPECT_EQ("simulation_test", input["simulation"]["name"].to<std::string>());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
