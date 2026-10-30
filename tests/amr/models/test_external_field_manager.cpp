#include "amr/physical_models/external_field_manager.hpp"
#include "core/models/external_field_updater_defs.hpp"

#include "initializer/data_provider.hpp"

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "gtest/gtest.h"

using namespace PHARE;
using namespace PHARE::core;


namespace
{
constexpr SimOpts mhd_opts{.dimension           = 2,
                           .interp_order        = 1,
                           .reconstruction_type = MHDOpts::ReconstructionType::Constant,
                           .slope_limiter_type  = MHDOpts::SlopeLimiterType::None,
                           .riemann_solver_type = MHDOpts::RiemannSolverType::Rusanov};

using MHDTypes  = PHARE_Types<mhd_opts>::MHD;
using Manager_t = solver::ExternalFieldManager<MHDTypes::VecField_t, MHDTypes::GridLayout_t>;

initializer::PHAREDict dictWithADipole()
{
    initializer::PHAREDict dict;
    auto& externalField            = dict["external_field"];
    externalField["type"]          = static_cast<int>(ExternalFieldUpdaterType::Dipole);
    externalField["position"]["x"] = -0.5;
    externalField["position"]["y"] = 0.5;
    externalField["moment"]["x"]   = 0.3;
    externalField["moment"]["y"]   = -0.2;
    externalField["radius"]        = 0.;
    return dict;
}

} // namespace


TEST(ExternalFieldManager, isDisabledWithoutExternalFieldKey)
{
    Manager_t manager{"model", initializer::PHAREDict{}};

    EXPECT_FALSE(manager.isEnabled());
    EXPECT_TRUE(manager.field().getRunTimeResourcesViewList().empty());
}


TEST(ExternalFieldManager, isEnabledByAnExternalFieldDict)
{
    Manager_t manager{"model", dictWithADipole()};

    EXPECT_TRUE(manager.isEnabled());
    EXPECT_EQ(manager.field().getRunTimeResourcesViewList().size(), 2u);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
