#include <SAMRAI/hier/BoxContainer.h>
#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/pdat/NodeGeometry.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>


#include "amr/data/field/field_geometry.hpp"
#include "amr/data/field/field_variable.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;

// this test checks that the FieldGeometry behaves the same
// as SAMRAI specific geometries for cases they are equivalent.
// namely, assuming a yee layout they chek that FieldGeometry
// associated with Ex is equivalent to the samrai cell geometry, and
// the the one associated with Ey is equivalent to the Node geometry.
// Note this is only done in 1D and cannot be done in 2 or 3D since there
// is no geometry in SAMRAI that can cope with the Yee layout in 2D or 3D.



using Field1D = Field<NdArrayVector<1>, HybridQuantity::Scalar>;

template<typename GridLayoutT, typename FieldImpl>
struct FieldGeometryParam
{
    FieldGeometryParam(std::string const& name, HybridQuantity::Scalar quantity,
                       SAMRAI::hier::Patch& patch_0, SAMRAI::hier::Patch& patch_1)
        : destinationFieldVariable{name + std::string("_0"), quantity}
        , sourceFieldVariable{name + std::string("_1"), quantity}
        , destinationFieldFactory{destinationFieldVariable.getPatchDataFactory()}
        , sourceFieldFactory{sourceFieldVariable.getPatchDataFactory()}
        , destinationPatch{patch_0}
        , sourcePatch{patch_1}
        , destinationFieldGeometry{destinationFieldFactory->getBoxGeometry(
              destinationPatch.getBox())}
        , sourceFieldGeometry{sourceFieldFactory->getBoxGeometry(sourcePatch.getBox())}
        , destinationFieldData{std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldImpl>>(
              destinationFieldFactory->allocate(destinationPatch))}
        , sourceFieldData{std::dynamic_pointer_cast<FieldData<GridLayoutT, FieldImpl>>(
              sourceFieldFactory->allocate(sourcePatch))}
    {
    }
    FieldVariable<GridLayoutT, FieldImpl> destinationFieldVariable;
    FieldVariable<GridLayoutT, FieldImpl> sourceFieldVariable;
    std::shared_ptr<SAMRAI::hier::PatchDataFactory> destinationFieldFactory;
    std::shared_ptr<SAMRAI::hier::PatchDataFactory> sourceFieldFactory;

    SAMRAI::hier::Patch& destinationPatch;
    SAMRAI::hier::Patch& sourcePatch;

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destinationFieldGeometry;
    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceFieldGeometry;

    std::shared_ptr<FieldData<GridLayoutT, FieldImpl>> destinationFieldData;
    std::shared_ptr<FieldData<GridLayoutT, FieldImpl>> sourceFieldData;
};




struct Patches1D
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};

    SAMRAI::hier::Box destinationDomainBox{SAMRAI::hier::Index(dim, 0),
                                           SAMRAI::hier::Index(dim, 10), blockId};
    SAMRAI::hier::Box sourceDomainBox{SAMRAI::hier::Index(dim, 5), SAMRAI::hier::Index(dim, 20),
                                      blockId};

    std::shared_ptr<SAMRAI::hier::PatchDescriptor> patchDescriptor{
        std::make_shared<SAMRAI::hier::PatchDescriptor>()};

    double dx                    = 0.01;
    double destinationPatchLower = 0.;
    double destinationPatchUpper = 0.1;

    double sourcePatchLower = 0.05;
    double sourcePatchUpper = 0.2;

    // create the patch geometry saying the patch does NOT touch a boundary
    // in none of the two directions of dim (==1)
    SAMRAI::hier::PatchGeometry::TwoDimBool touchesRegular{dim, false};

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> destinationPatchGeometry{
        std::make_shared<SAMRAI::geom::CartesianPatchGeometry>(
            SAMRAI::hier::IntVector::getOne(dim), touchesRegular, blockId, &dx,
            &destinationPatchLower, &destinationPatchUpper)};

    std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> sourcePatchGeometry{
        std::make_shared<SAMRAI::geom::CartesianPatchGeometry>(
            SAMRAI::hier::IntVector::getOne(dim), touchesRegular, blockId, &dx, &sourcePatchLower,
            &sourcePatchUpper)};

    SAMRAI::hier::Patch destinationPatch{destinationDomainBox, patchDescriptor};
    SAMRAI::hier::Patch sourcePatch{sourceDomainBox, patchDescriptor};

    Patches1D()
    {
        destinationPatch.setPatchGeometry(destinationPatchGeometry);
        sourcePatch.setPatchGeometry(sourcePatchGeometry);
    }
};




template<typename T>
struct FieldGeometry1D : public ::testing::Test
{
};




TYPED_TEST_SUITE_P(FieldGeometry1D);

TYPED_TEST_P(FieldGeometry1D, IsSameAsCellGeometryForEx)
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};
    Patches1D patches1d;

    auto& destinationPatch = patches1d.destinationPatch;
    auto& sourcePatch      = patches1d.sourcePatch;

    TypeParam param{"Ex", HybridQuantity::Scalar::Ex, destinationPatch, sourcePatch};

    auto& destinationLayout = param.destinationFieldData->gridLayout;
    auto centering          = destinationLayout.centering(HybridQuantity::Scalar::Ex);
    auto ghosts             = SAMRAI::hier::IntVector::getZero(dim);
    ghosts[0]               = destinationLayout.nbrGhosts(centering[0]);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destinationCellGeometry
        = std::make_shared<SAMRAI::pdat::CellGeometry>(destinationPatch.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceCellGeometry
        = std::make_shared<SAMRAI::pdat::CellGeometry>(sourcePatch.getBox(), ghosts);

    int lower = 6;
    int upper = 11;

    // We want here to get the ghost for interp order >=2
    if (destinationLayout.interp_order >= 2)
    {
        upper = 12;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};
    SAMRAI::hier::Box fillBox{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};

    SAMRAI::hier::Box restrictBox1{SAMRAI::hier::Index{dim, 7}, SAMRAI::hier::Index{dim, 8},
                                   blockId};
    SAMRAI::hier::Box restrictBox2{SAMRAI::hier::Index{dim, 10}, SAMRAI::hier::Index{dim, 11},
                                   blockId};

    std::array<bool, 2> overwritePossibility{{false, true}};

    SAMRAI::hier::BoxContainer noRestrictBoxes{};

    SAMRAI::hier::BoxContainer someRestrictBoxes{restrictBox1};
    someRestrictBoxes.push_back(restrictBox2);

    std::array<SAMRAI::hier::BoxContainer, 2> restrictBoxesList{
        {noRestrictBoxes, someRestrictBoxes}};

    for (auto const& restrictBoxes : restrictBoxesList)
    {
        for (auto overwriteInterior : overwritePossibility)
        {
            SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(dim)};

            auto fieldOverlap = std::dynamic_pointer_cast<FieldOverlap>(
                param.destinationFieldGeometry->calculateOverlap(
                    *param.sourceFieldGeometry, srcMask, fillBox, overwriteInterior, transformation,
                    restrictBoxes));

            auto cellOverlap = std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(
                destinationCellGeometry->calculateOverlap(*sourceCellGeometry, srcMask, fillBox,
                                                          overwriteInterior, transformation,
                                                          restrictBoxes));

            EXPECT_THAT(fieldOverlap->isOverlapEmpty(), Eq(cellOverlap->isOverlapEmpty()));

            auto fieldDestBoxes = fieldOverlap->getDestinationBoxContainer();
            auto cellDestBoxes  = cellOverlap->getDestinationBoxContainer();

            EXPECT_THAT(fieldDestBoxes, Eq(cellDestBoxes));
        }
    }
}




TYPED_TEST_P(FieldGeometry1D, IsSameAsNodeGeometryForEy)
{
    SAMRAI::tbox::Dimension dim{1};
    SAMRAI::hier::BlockId blockId{0};

    Patches1D patch1d;

    auto& destinationPatch = patch1d.destinationPatch;
    auto& sourcePatch      = patch1d.sourcePatch;

    TypeParam param{"Ey", HybridQuantity::Scalar::Ey, destinationPatch, sourcePatch};

    auto& destinationLayout = param.destinationFieldData->gridLayout;

    auto centering = destinationLayout.centering(HybridQuantity::Scalar::Ey);

    auto ghosts = SAMRAI::hier::IntVector::getZero(dim);

    ghosts[0] = destinationLayout.nbrGhosts(centering[0]);


    std::shared_ptr<SAMRAI::hier::BoxGeometry> node0Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(destinationPatch.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> node1Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(sourcePatch.getBox(), ghosts);

    int lower = 6;
    int upper = 11;

    // We want here to get the ghost for interp order >=2
    if (destinationLayout.interp_order >= 2)
    {
        upper = 12;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};
    SAMRAI::hier::Box fillBox{SAMRAI::hier::Index{dim, lower}, SAMRAI::hier::Index{dim, upper},
                              blockId};

    SAMRAI::hier::Box restrictBox1{SAMRAI::hier::Index{dim, 7}, SAMRAI::hier::Index{dim, 8},
                                   blockId};
    SAMRAI::hier::Box restrictBox2{SAMRAI::hier::Index{dim, 10}, SAMRAI::hier::Index{dim, 11},
                                   blockId};


    SAMRAI::hier::BoxContainer noRestrictBoxes{};

    SAMRAI::hier::BoxContainer someRestrictBoxes{restrictBox1};
    someRestrictBoxes.push_back(restrictBox2);

    std::array<SAMRAI::hier::BoxContainer, 2> restrictBoxesList{
        {noRestrictBoxes, someRestrictBoxes}};

    std::array<bool, 2> overwritePossibility{{false, true}};

    for (auto const& restrictBoxes : restrictBoxesList)
    {
        for (auto overwriteInterior : overwritePossibility)
        {
            SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(dim)};

            auto fieldOverlap = std::dynamic_pointer_cast<FieldOverlap>(
                param.destinationFieldGeometry->calculateOverlap(
                    *param.sourceFieldGeometry, srcMask, fillBox, overwriteInterior, transformation,
                    restrictBoxes));

            auto nodeOverlap = std::dynamic_pointer_cast<SAMRAI::pdat::NodeOverlap>(
                node0Geom->calculateOverlap(*node1Geom, srcMask, fillBox, overwriteInterior,
                                            transformation, restrictBoxes));

            EXPECT_THAT(fieldOverlap->isOverlapEmpty(), Eq(nodeOverlap->isOverlapEmpty()));

            auto fieldDestBoxes = fieldOverlap->getDestinationBoxContainer();
            auto nodeDestBoxes  = nodeOverlap->getDestinationBoxContainer();

            EXPECT_THAT(fieldDestBoxes, Eq(nodeDestBoxes));
        }
    }
}



REGISTER_TYPED_TEST_SUITE_P(FieldGeometry1D, IsSameAsCellGeometryForEx, IsSameAsNodeGeometryForEy);

using FieldGeometryTest1DOrder1 = FieldGeometryParam<GridLayout<GridLayoutImplYee<1, 1>>, Field1D>;
using FieldGeometryTest1DOrder2 = FieldGeometryParam<GridLayout<GridLayoutImplYee<1, 2>>, Field1D>;
using FieldGeometryTest1DOrder3 = FieldGeometryParam<GridLayout<GridLayoutImplYee<1, 3>>, Field1D>;

using FieldGeometry1DTestList
    = ::testing::Types<FieldGeometryTest1DOrder1, FieldGeometryTest1DOrder2,
                       FieldGeometryTest1DOrder3>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestWithOrderFrom1To3That, FieldGeometry1D, FieldGeometry1DTestList);


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
