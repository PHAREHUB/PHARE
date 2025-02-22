
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include "amr/data/field/field_variable.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

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



template<std::size_t dim, std::size_t interporder>
using FieldVar = FieldVariable<GridLayout<GridLayoutImplYee<dim, interporder>>,
                               Grid<NdArrayVector<dim, floater_t<4>>, HybridQuantity::Scalar>>;


// The interp order is of no importance to know if a quantity
// lives on or inside a patch boundary



using TestWithQuantityThatLivesOnPatchBoundary1D = FieldVariableTest;
using TestWithQuantityThatLivesOnPatchBoundary2D = FieldVariableTest;
using TestWithQuantityThatLivesOnPatchBoundary3D = FieldVariableTest;




TEST_P(TestWithQuantityThatLivesOnPatchBoundary1D, ThatActualDataLivesOnPatchBoundary)
{
    auto constexpr dim    = 1;
    auto constexpr interp = 1;

    auto fieldVariable = std::make_shared<FieldVar<dim, interp>>(param.qtyName, param.qty);
    EXPECT_TRUE(fieldVariable->dataLivesOnPatchBorder());
}


TEST_P(TestWithQuantityThatLivesOnPatchBoundary2D, ThatActualDataLivesOnPatchBoundary)
{
    auto constexpr dim    = 2;
    auto constexpr interp = 1;

    auto fieldVariable = std::make_shared<FieldVar<dim, interp>>(param.qtyName, param.qty);
    EXPECT_TRUE(fieldVariable->dataLivesOnPatchBorder());
}


TEST_P(TestWithQuantityThatLivesOnPatchBoundary3D, ThatActualDataLivesOnPatchBoundary)
{
    auto constexpr dim    = 3;
    auto constexpr interp = 1;

    auto fieldVariable = std::make_shared<FieldVar<dim, interp>>(param.qtyName, param.qty);
    EXPECT_TRUE(fieldVariable->dataLivesOnPatchBorder());
}



using TestWithQuantityThatLivesInsidePatchBoundary1D = FieldVariableTest;
using TestWithQuantityThatLivesInsidePatchBoundary2D = FieldVariableTest;
using TestWithQuantityThatLivesInsidePatchBoundary3D = FieldVariableTest;


TEST_P(TestWithQuantityThatLivesInsidePatchBoundary1D, ThatActualDataLivesInsidePatchBoundary)
{
    auto constexpr dim    = 1;
    auto constexpr interp = 1;

    auto fieldVariable = std::make_shared<FieldVar<dim, interp>>(param.qtyName, param.qty);
    EXPECT_FALSE(fieldVariable->dataLivesOnPatchBorder());
}


TEST_P(TestWithQuantityThatLivesInsidePatchBoundary2D, ThatActualDataLivesInsidePatchBoundary)
{
    auto constexpr dim    = 2;
    auto constexpr interp = 1;

    auto fieldVariable = std::make_shared<FieldVar<dim, interp>>(param.qtyName, param.qty);
    EXPECT_FALSE(fieldVariable->dataLivesOnPatchBorder());
}


TEST_P(TestWithQuantityThatLivesInsidePatchBoundary3D, ThatActualDataLivesInsidePatchBoundary)
{
    auto constexpr dim    = 3;
    auto constexpr interp = 1;

    auto fieldVariable = std::make_shared<FieldVar<dim, interp>>(param.qtyName, param.qty);
    EXPECT_FALSE(fieldVariable->dataLivesOnPatchBorder());
}


// the definition of which variable lives on or inside patch boundaries
// is hard-coded for the Yee Layout Implementation @ 1D
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


// the same for 2D
std::map<std::string, HybridQuantity::Scalar> On2DPatchBoundaryQties
    = {{"Bx", HybridQuantity::Scalar::Bx},   {"By", HybridQuantity::Scalar::By},
       {"Ex", HybridQuantity::Scalar::Ex},   {"Ey", HybridQuantity::Scalar::Ey},
       {"Ez", HybridQuantity::Scalar::Ez},   {"Jx", HybridQuantity::Scalar::Jx},
       {"Jy", HybridQuantity::Scalar::Jy},   {"Jz", HybridQuantity::Scalar::Jz},
       {"rho", HybridQuantity::Scalar::rho}, {"Vx", HybridQuantity::Scalar::Vx},
       {"Vy", HybridQuantity::Scalar::Vy},   {"Vz", HybridQuantity::Scalar::Vz},
       {"P", HybridQuantity::Scalar::P}};

std::map<std::string, HybridQuantity::Scalar> Inside2DPatchBoundaryQties
    = {{"Bz", HybridQuantity::Scalar::Bz}};


// and at 3D, everybody is on the patch boundary
std::map<std::string, HybridQuantity::Scalar> On3DPatchBoundaryQties
    = {{"Bx", HybridQuantity::Scalar::Bx}, {"By", HybridQuantity::Scalar::By},
       {"Bz", HybridQuantity::Scalar::Bz}, {"Ex", HybridQuantity::Scalar::Ex},
       {"Ey", HybridQuantity::Scalar::Ey}, {"Ez", HybridQuantity::Scalar::Ez},
       {"Jx", HybridQuantity::Scalar::Jx}, {"Jy", HybridQuantity::Scalar::Jy},
       {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
       {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
       {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};

std::map<std::string, HybridQuantity::Scalar> Inside3DPatchBoundaryQties = {};

// because "Inside3DPatchBoundaryQties" is empty so the associated test is not instantiated
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(TestWithQuantityThatLivesInsidePatchBoundary3D);



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

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesOnPatchBoundary2D,
                         ::testing::ValuesIn(createParams(On2DPatchBoundaryQties)));

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesOnPatchBoundary3D,
                         ::testing::ValuesIn(createParams(On3DPatchBoundaryQties)));


INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesInsidePatchBoundary1D,
                         ::testing::ValuesIn(createParams(Inside1DPatchBoundaryQties)));

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesInsidePatchBoundary2D,
                         ::testing::ValuesIn(createParams(Inside2DPatchBoundaryQties)));

INSTANTIATE_TEST_SUITE_P(FieldVariable, TestWithQuantityThatLivesInsidePatchBoundary3D,
                         ::testing::ValuesIn(createParams(Inside3DPatchBoundaryQties)));



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
