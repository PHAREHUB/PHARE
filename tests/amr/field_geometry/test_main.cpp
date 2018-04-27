#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/pdat/NodeGeometry.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/field/field_geometry.h"
#include "data/field/field_variable.h"

using testing::Eq;

namespace PHARE
{
using Field1D = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

template<Layout layout, std::size_t dim, std::size_t interpOrder, typename FieldImpl>
struct FieldGeometryParam
{
    FieldGeometryParam(std::string const& name, HybridQuantity::Scalar quantity,
                       SAMRAI::hier::Patch& patch_0, SAMRAI::hier::Patch& patch_1)
        : field0Variable{name + std::string("_0"), true, quantity}
        , field1Variable{name + std::string("_1"), true, quantity}
        , field0Factory{field0Variable.getPatchDataFactory()}
        , field1Factory{field1Variable.getPatchDataFactory()}
        , patch0{patch_0}
        , patch1{patch_1}
        , field0Geom{field0Factory->getBoxGeometry(patch0.getBox())}
        , field1Geom{field1Factory->getBoxGeometry(patch1.getBox())}
        , field0Data{std::dynamic_pointer_cast<FieldData<layout, dim, interpOrder, FieldImpl>>(
              field0Factory->allocate(patch0))}
        , field1Data{std::dynamic_pointer_cast<FieldData<layout, dim, interpOrder, FieldImpl>>(
              field1Factory->allocate(patch1))}
    {
    }
    FieldVariable<layout, dim, interpOrder, FieldImpl> field0Variable;
    FieldVariable<layout, dim, interpOrder, FieldImpl> field1Variable;
    std::shared_ptr<SAMRAI::hier::PatchDataFactory> field0Factory;
    std::shared_ptr<SAMRAI::hier::PatchDataFactory> field1Factory;

    SAMRAI::hier::Patch& patch0;
    SAMRAI::hier::Patch& patch1;

    std::shared_ptr<SAMRAI::hier::BoxGeometry> field0Geom;
    std::shared_ptr<SAMRAI::hier::BoxGeometry> field1Geom;

    std::shared_ptr<FieldData<layout, dim, interpOrder, FieldImpl>> field0Data;
    std::shared_ptr<FieldData<layout, dim, interpOrder, FieldImpl>> field1Data;
};

struct PatchGeometry1D
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};


    SAMRAI::hier::Box box0{SAMRAI::hier::Index(dim, 0), SAMRAI::hier::Index(dim, 10), blockId};
    SAMRAI::hier::Box box1{SAMRAI::hier::Index(dim, 5), SAMRAI::hier::Index(dim, 20), blockId};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    double dx{0.01};
    double patch0_lo{0.};
    double patch0_hi{0.1};

    double patch1_lo{0.1};
    double patch1_hi{0.2};


    SAMRAI::hier::PatchGeometry::TwoDimBool touchesRegular{dim, false};

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch0Geom{
        std::make_shared<SAMRAI::geom::CartesianPatchGeometry>(SAMRAI::hier::IntVector::getOne(dim),
                                                               touchesRegular, blockId, &dx,
                                                               &patch0_lo, &patch0_hi)};

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> patch1Geom{
        std::make_shared<SAMRAI::geom::CartesianPatchGeometry>(SAMRAI::hier::IntVector::getOne(dim),
                                                               touchesRegular, blockId, &dx,
                                                               &patch1_lo, &patch1_hi)};



    SAMRAI::hier::Patch patch0{box0, patchDescriptor};
    SAMRAI::hier::Patch patch1{box1, patchDescriptor};

    PatchGeometry1D()
    {
        patch0.setPatchGeometry(patch0Geom);
        patch1.setPatchGeometry(patch1Geom);
    }
};

template<typename T>
struct FieldGeometry1D : public ::testing::Test
{
};

TYPED_TEST_CASE_P(FieldGeometry1D);

TYPED_TEST_P(FieldGeometry1D, IsSameAsCellGeometryForEx)
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};


    PatchGeometry1D patch1d;

    auto& patch0 = patch1d.patch0;
    auto& patch1 = patch1d.patch1;

    TypeParam param{"Ex", HybridQuantity::Scalar::Ex, patch0, patch1};

    auto& layout0 = param.field0Data->gridLayout;

    auto centering = layout0.centering(HybridQuantity::Scalar::Ex);


    auto ghosts = SAMRAI::hier::IntVector::getZero(dim);

    ghosts[0] = layout0.nbrGhostNodes(centering[0]);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> cell0Geom
        = std::make_shared<SAMRAI::pdat::CellGeometry>(patch0.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> cell1Geom
        = std::make_shared<SAMRAI::pdat::CellGeometry>(patch1.getBox(), ghosts);

    int lower = 6;
    int upper = 9;

    if (layout0.order() == 2)
    {
        upper = 15;
    }
    else if (layout0.order() == 3)
    {
        upper = 20;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};
    SAMRAI::hier::Box fillBox{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};

    bool overwriteInterior{true};

    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(dim)};

    auto fieldOverlap
        = std::dynamic_pointer_cast<FieldOverlap<1>>(param.field0Geom->calculateOverlap(
            *param.field1Geom, srcMask, fillBox, overwriteInterior, transformation));

    auto cellOverlap
        = std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(cell0Geom->calculateOverlap(
            *cell1Geom, srcMask, fillBox, overwriteInterior, transformation));

    EXPECT_THAT(fieldOverlap->isOverlapEmpty(), Eq(cellOverlap->isOverlapEmpty()));

    auto fieldDestBoxes = fieldOverlap->getDestinationBoxContainer();
    auto cellDestBoxes  = cellOverlap->getDestinationBoxContainer();

    EXPECT_THAT(fieldDestBoxes, Eq(cellDestBoxes));
}

TYPED_TEST_P(FieldGeometry1D, IsSameAsNodeGeometryForEy)
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};


    PatchGeometry1D patch1d;

    auto& patch0 = patch1d.patch0;
    auto& patch1 = patch1d.patch1;

    FieldGeometryParam<Layout::Yee, 1, 1, Field1D> param{"Ey", HybridQuantity::Scalar::Ey, patch0,
                                                         patch1};

    auto& layout0 = param.field0Data->gridLayout;

    auto centering = layout0.centering(HybridQuantity::Scalar::Ey);

    auto ghosts = SAMRAI::hier::IntVector::getZero(dim);

    ghosts[0] = layout0.nbrGhostNodes(centering[0]);


    std::shared_ptr<SAMRAI::hier::BoxGeometry> node0Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(patch0.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> node1Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(patch1.getBox(), ghosts);

    int lower = 6;
    int upper = 9;

    if (layout0.order() == 2)
    {
        upper = 15;
    }
    else if (layout0.order() == 3)
    {
        upper = 20;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};
    SAMRAI::hier::Box fillBox{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};

    bool overwriteInterior{true};

    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(dim)};

    auto fieldOverlap
        = std::dynamic_pointer_cast<FieldOverlap<1>>(param.field0Geom->calculateOverlap(
            *param.field1Geom, srcMask, fillBox, overwriteInterior, transformation));

    auto nodeOverlap
        = std::dynamic_pointer_cast<SAMRAI::pdat::NodeOverlap>(node0Geom->calculateOverlap(
            *node1Geom, srcMask, fillBox, overwriteInterior, transformation));

    EXPECT_THAT(fieldOverlap->isOverlapEmpty(), Eq(nodeOverlap->isOverlapEmpty()));

    auto fieldDestBoxes = fieldOverlap->getDestinationBoxContainer();
    auto nodeDestBoxes  = nodeOverlap->getDestinationBoxContainer();

    EXPECT_THAT(fieldDestBoxes, Eq(nodeDestBoxes));
}
REGISTER_TYPED_TEST_CASE_P(FieldGeometry1D, IsSameAsCellGeometryForEx, IsSameAsNodeGeometryForEy);

using FieldGeometryTest1DOrder1 = FieldGeometryParam<Layout::Yee, 1, 1, Field1D>;
using FieldGeometryTest1DOrder2 = FieldGeometryParam<Layout::Yee, 1, 2, Field1D>;
using FieldGeometryTest1DOrder3 = FieldGeometryParam<Layout::Yee, 1, 3, Field1D>;

using FieldGeometry1DTestList
    = ::testing::Types<FieldGeometryTest1DOrder1, FieldGeometryTest1DOrder2,
                       FieldGeometryTest1DOrder3>;

INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, FieldGeometry1D, FieldGeometry1DTestList);

} // namespace PHARE

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
