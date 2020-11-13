#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include "core/data/field/field.h"
#include "amr/data/field/field_variable.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/ndarray/ndarray_vector.h"

#include <string>
#include <map>


using namespace PHARE::core;
using namespace PHARE::amr;

struct FieldVariableTestParam
{
public:
    FieldVariableTestParam() = default;
    FieldVariableTestParam(std::string const& name, HybridQuantity::Scalar quantity)
        : qtyName{name}
        , qty{quantity}
    {
    }
    std::string qtyName;
    HybridQuantity::Scalar qty;
};

struct FieldVariableTest : public ::testing::TestWithParam<FieldVariableTestParam>
{
    FieldVariableTest() = default;
    void SetUp() override { param = GetParam(); }
    virtual ~FieldVariableTest() = default;

    FieldVariableTestParam param;
};


using TestWithQuantityThatLivesOnPatchBoundary1D = FieldVariableTest;

template <std::size_t dim, std::size_t interporder>
using FV = FieldVariable<GridLayout<GridLayoutImplYee<dim, interporder>>,
                         Field<NdArrayVector<dim>, HybridQuantity::Scalar>>;

TEST_P(TestWithQuantityThatLivesOnPatchBoundary1D, ThatActualDataLivesOnPatchBoundary)
{
    auto fieldVariable  = std::make_shared<FV<1,1>>(
            param.qtyName, param.qty);

    EXPECT_TRUE(fieldVariable->dataLivesOnPatchBorder());
}

using TestWithQuantityThatLivesInsidePatchBoundary1D = FieldVariableTest;

TEST_P(TestWithQuantityThatLivesInsidePatchBoundary1D, ThatActualDataLivesInsidePatchBoundary)
{
    auto fieldVariable
        = std::make_shared<FV<1,1>>(
            param.qtyName, param.qty);

    EXPECT_FALSE(fieldVariable->dataLivesOnPatchBorder());
}


// the definition of which variable lives on or inside patch boundaries
// depends on the dimension and here is hard-coded for the Yee Layout Implementation
std::map<std::string, HybridQuantity::Scalar> On1DPatchBoundaryQties
    = {{"Bx", HybridQuantity::Scalar::Bx}, {"Ey", HybridQuantity::Scalar::Ey},
       {"Ez", HybridQuantity::Scalar::Ez}, {"Jy", HybridQuantity::Scalar::Jy},
       {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
       {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
       {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};

std::map<std::string, HybridQuantity::Scalar> Inside1DPatchBoundaryQties
    = {{"By", HybridQuantity::Scalar::By},
       {"Bz", HybridQuantity::Scalar::Bz},
       {"Ex", HybridQuantity::Scalar::Ex},
       {"Jx", HybridQuantity::Scalar::Jx}};


std::vector<FieldVariableTestParam>
createParams(std::map<std::string, HybridQuantity::Scalar> const& qtyMap)
{
    std::vector<FieldVariableTestParam> params;
    for (auto const& qtyPair : qtyMap)
    {
        params.emplace_back(qtyPair.first, qtyPair.second);
    }
    return params;
}

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesOnPatchBoundary1D,
                         ::testing::ValuesIn(createParams(On1DPatchBoundaryQties)));

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesInsidePatchBoundary1D,
                         ::testing::ValuesIn(createParams(Inside1DPatchBoundaryQties)));


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
