#include <type_traits>

#include "core/data/field/field.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/models/mhd_state.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/vecfield/vecfield.hpp"

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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
