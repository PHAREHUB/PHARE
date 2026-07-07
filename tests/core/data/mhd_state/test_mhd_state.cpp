#include "core/data/field/field.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"
#include "core/numerics/primite_conservative_converter/mhd_conversion.hpp"

#include "initializer/data_provider.hpp"

#include "gtest/gtest.h"

#include "tests/core/data/mhd_state/init_functions.hpp"

using namespace PHARE::core;
using namespace PHARE::initializer;
using namespace PHARE::initializer::test_fn::func_1d;

using Field_t    = Field<1, MHDQuantity::Scalar>;
using VecField1D = VecField<Field_t, MHDQuantity>;


PHAREDict getDict()
{
    using initfunc = InitFunction<1>;
    PHAREDict dict;

    dict["name"] = std::string("state");

    dict["density"]["initializer"] = static_cast<initfunc>(density);

    dict["velocity"]["initializer"]["x_component"] = static_cast<initfunc>(vx);
    dict["velocity"]["initializer"]["y_component"] = static_cast<initfunc>(vy);
    dict["velocity"]["initializer"]["z_component"] = static_cast<initfunc>(vz);

    dict["magnetic"]["initializer"]["x_component"] = static_cast<initfunc>(bx);
    dict["magnetic"]["initializer"]["y_component"] = static_cast<initfunc>(by);
    dict["magnetic"]["initializer"]["z_component"] = static_cast<initfunc>(bz);

    dict["pressure"]["initializer"] = static_cast<initfunc>(pressure);

    dict["to_conservative_init"]["heat_capacity_ratio"] = 5. / 3.;

    return dict;
}

struct AnMHDState : public ::testing::Test
{
    MHDState<VecField1D> state{getDict()};
    virtual ~AnMHDState();
};

AnMHDState::~AnMHDState() {}

TEST_F(AnMHDState, noUsableFieldsUponConstruction)
{
    EXPECT_FALSE(state.isUsable());
}

TEST_F(AnMHDState, fieldsAreSettable)
{
    EXPECT_TRUE(state.isSettable());
}

TEST_F(AnMHDState, hasTupleResourceList)
{
    auto resources              = state.getCompileTimeResourcesViewList();
    [[maybe_unused]] auto& rho  = std::get<0>(resources);
    [[maybe_unused]] auto& v    = std::get<1>(resources);
    [[maybe_unused]] auto& b    = std::get<2>(resources);
    [[maybe_unused]] auto& p    = std::get<3>(resources);
    [[maybe_unused]] auto& rhov = std::get<4>(resources);
    [[maybe_unused]] auto& etot = std::get<5>(resources);
    [[maybe_unused]] auto& j    = std::get<6>(resources);
    [[maybe_unused]] auto& e    = std::get<7>(resources);

    EXPECT_FALSE(rho.isUsable());
    EXPECT_FALSE(v.isUsable());
    EXPECT_FALSE(b.isUsable());
    EXPECT_FALSE(p.isUsable());
    EXPECT_FALSE(rhov.isUsable());
    EXPECT_FALSE(etot.isUsable());
    EXPECT_FALSE(j.isUsable());
    EXPECT_FALSE(e.isUsable());
}

TEST(MHDConversion, reconstructsTotalMagneticField)
{
    auto const& [bx, by, bz] = totalMagneticComponents(1.0, -2.0, 0.5, 0.25, 1.5, -1.0);

    EXPECT_DOUBLE_EQ(bx, 1.25);
    EXPECT_DOUBLE_EQ(by, -0.5);
    EXPECT_DOUBLE_EQ(bz, -0.5);
}

TEST(MHDConversion, convertsReducedAndTotalEnergy)
{
    auto constexpr gamma = 5. / 3.;
    auto constexpr rho   = 4.0;
    auto constexpr vx    = 1.0;
    auto constexpr vy    = -0.5;
    auto constexpr vz    = 0.25;
    auto constexpr b1x   = 0.3;
    auto constexpr b1y   = -0.4;
    auto constexpr b1z   = 0.5;
    auto constexpr b0x   = 1.5;
    auto constexpr b0y   = -0.25;
    auto constexpr b0z   = 0.75;
    auto constexpr p     = 2.0;

    auto reducedEnergy     = eosPToEtot1(gamma, rho, vx, vy, vz, b1x, b1y, b1z, p);
    auto totalEnergy       = etot1ToEtot(reducedEnergy, b1x, b1y, b1z, b0x, b0y, b0z);
    auto reducedEnergyBack = etotToEtot1(totalEnergy, b1x, b1y, b1z, b0x, b0y, b0z);
    auto pressureBack      = eosEtot1ToP(gamma, rho, vx, vy, vz, b1x, b1y, b1z, reducedEnergyBack);
    auto const& [bx, by, bz] = totalMagneticComponents(b1x, b1y, b1z, b0x, b0y, b0z);
    auto totalEnergyDirect   = eosPToEtot(gamma, rho, vx, vy, vz, bx, by, bz, p);

    EXPECT_DOUBLE_EQ(reducedEnergyBack, reducedEnergy);
    EXPECT_DOUBLE_EQ(pressureBack, p);
    EXPECT_DOUBLE_EQ(totalEnergy, totalEnergyDirect);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
