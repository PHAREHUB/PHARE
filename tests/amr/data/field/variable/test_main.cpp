#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include "core/data/field/field.h"
#include "amr/data/field/field_variable.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/ndarray/ndarray_vector.h"



using namespace PHARE::core;
using namespace PHARE::amr;

struct FieldVariableTestParam1D
{
public:
    FieldVariableTestParam1D() = default;
    FieldVariableTestParam1D(std::string const& name, HybridQuantity::Scalar quantity)
        : qtyName{name}
        , qty{quantity}
    {
    }

    std::string qtyName;
    HybridQuantity::Scalar qty;
};

struct FieldVariableTest1D : public ::testing::TestWithParam<FieldVariableTestParam1D>
{
    FieldVariableTest1D() = default;
    void SetUp() override { param = GetParam(); }
    virtual ~FieldVariableTest1D() = default;

    FieldVariableTestParam1D param;
};


using TestWithQuantityThatLivesOnPatchBoundary1D = FieldVariableTest1D;

TEST_P(TestWithQuantityThatLivesOnPatchBoundary1D, ThatActualDataLivesOnPatchBoundary)
{
    auto fieldVariable
        = std::make_shared<FieldVariable<GridLayout<GridLayoutImplYee<1, 1>>,
                                         Field<NdArrayVector<1>, HybridQuantity::Scalar>>>(
            param.qtyName, param.qty);

    EXPECT_TRUE(fieldVariable->dataLivesOnPatchBorder());
}

using TestWithQuantityThatLivesInsidePatchBoundary1D = FieldVariableTest1D;

TEST_P(TestWithQuantityThatLivesInsidePatchBoundary1D, ThatActualDataLivesInsidePatchBoundary)
{
    auto fieldVariable
        = std::make_shared<FieldVariable<GridLayout<GridLayoutImplYee<1, 1>>,
                                         Field<NdArrayVector<1>, HybridQuantity::Scalar>>>(
            param.qtyName, param.qty);

    EXPECT_FALSE(fieldVariable->dataLivesOnPatchBorder());
}

std::map<std::string, HybridQuantity::Scalar> quantityThatLivesOnPatchBoundary1D
    = {{"Bx", HybridQuantity::Scalar::Bx}, {"Ey", HybridQuantity::Scalar::Ey},
       {"Ez", HybridQuantity::Scalar::Ez}, {"Jy", HybridQuantity::Scalar::Jy},
       {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
       {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
       {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};

std::map<std::string, HybridQuantity::Scalar> quantityThatLivesInsidePatchBoundary1D
    = {{"By", HybridQuantity::Scalar::By},
       {"Bz", HybridQuantity::Scalar::Bz},
       {"Ex", HybridQuantity::Scalar::Ex},
       {"Jx", HybridQuantity::Scalar::Jx}};

std::vector<FieldVariableTestParam1D>
createParams(std::map<std::string, HybridQuantity::Scalar> const& qtyMap)
{
    std::vector<FieldVariableTestParam1D> params;

    for (auto const& qtyPair : qtyMap)
    {
        params.emplace_back(qtyPair.first, qtyPair.second);
    }

    return params;
}

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesOnPatchBoundary1D,
                         ::testing::ValuesIn(createParams(quantityThatLivesOnPatchBoundary1D)));

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesInsidePatchBoundary1D,
                         ::testing::ValuesIn(createParams(quantityThatLivesInsidePatchBoundary1D)));


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
