// -------------------------------------
//   FieldLinearTimeInterpolate test
// -------------------------------------

#include <type_traits>

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "amr/data/field/time_interpolate/field_linear_time_interpolate.h"

#include "core/data/field/field.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/hybrid/hybrid_quantities.h"
#include "amr/resources_manager/amr_utils.h"



using namespace PHARE::core;
using namespace PHARE::amr;

std::size_t constexpr dim         = 1;
std::size_t constexpr interpOrder = 1;


using GridYee = GridLayout<GridLayoutImplYee<dim, interpOrder>>;

using Field1D = Field<NdArrayVector<1>, HybridQuantity::Scalar>;

using FieldDataT = FieldData<GridYee, Field1D>;


struct aFieldLinearTimeInterpolate : public ::testing::Test
{
    FieldLinearTimeInterpolate<GridYee, Field1D> timeOp{};

    HybridQuantity::Scalar qty{HybridQuantity::Scalar::Bx};

    SAMRAI::tbox::Dimension dimension{dim};

    SAMRAI::hier::BlockId block0{0};

    static int countLocal;

    SAMRAI::hier::Box domain{SAMRAI::hier::Index{dimension, 0}, SAMRAI::hier::Index{dimension, 100},
                             block0, SAMRAI::hier::LocalId{++countLocal},
                             SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank()};

    std::string const fieldName{"Bx"};

    SAMRAI::hier::IntVector ghost{dimension, 5};

    std::array<double, dim> dl{{0.01}};

    std::array<std::uint32_t, dim> nbrCells{{static_cast<std::uint32_t>(domain.numberCells(dirX))}};

    Point<double, dim> origin{0.};

    std::shared_ptr<FieldDataT> srcOld;
    std::shared_ptr<FieldDataT> srcNew;
    std::shared_ptr<FieldDataT> destNew;

    double srcFunc(double x, double t) { return x + t; }


    aFieldLinearTimeInterpolate()
        : srcOld{std::make_shared<FieldDataT>(domain, ghost, fieldName, dl, nbrCells, origin, qty)}
        , srcNew{std::make_shared<FieldDataT>(domain, ghost, fieldName, dl, nbrCells, origin, qty)}
        , destNew{std::make_shared<FieldDataT>(domain, ghost, fieldName, dl, nbrCells, origin, qty)}
    {
        // first we initialize time for each patch data
        // we test that alpha = 0, so that the interpolate time is oldTime

        double oldTime = 0.;
        double newTime = 0.5;

        srcOld->setTime(oldTime);
        srcNew->setTime(newTime);




        auto& layout = srcOld->gridLayout;

        auto& srcFieldOld = srcOld->field;
        auto& srcFieldNew = srcNew->field;


        auto iStart = layout.ghostStartIndex(qty, Direction::X);
        auto iEnd   = layout.ghostEndIndex(qty, Direction::X);

        for (auto ix = iStart; ix <= iEnd; ++ix)
        {
            auto position = layout.fieldNodeCoordinates(srcFieldOld, origin, ix);

            srcFieldOld(ix) = srcFunc(position[dirX], oldTime);

            srcFieldNew(ix) = srcFunc(position[dirX], newTime);
        }
    }


    auto zeroTransformation()
    {
        return SAMRAI::hier::Transformation(SAMRAI::hier::Transformation::NO_ROTATE,
                                            SAMRAI::hier::IntVector::getZero(dimension),
                                            SAMRAI::hier::BlockId(0), SAMRAI::hier::BlockId(0));
    }
};

int aFieldLinearTimeInterpolate::countLocal = 0;

TEST_F(aFieldLinearTimeInterpolate, giveOldSrcForAlphaZero)
{
    double interpolateTime = 0.;
    destNew->setTime(interpolateTime);

    auto& layout = srcOld->gridLayout;


    auto& srcFieldOld = srcOld->field;

    auto& destField = destNew->field;

    auto zero_transformation{zeroTransformation()};
    SAMRAI::hier::BoxContainer ghost_cntnr;
    FieldOverlap<dim> overlap{ghost_cntnr, zero_transformation};

    timeOp.timeInterpolate(*destNew, domain, overlap, *srcOld, *srcNew);


    bool const withGhost{true};
    auto box = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(domain, qty, layout,
                                                                          !withGhost);

    auto ghostBox = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(domain, qty, layout,
                                                                               withGhost);

    auto localBox = AMRToLocal(static_cast<std::add_const_t<decltype(box)>>(box), ghostBox);

    auto iStart = localBox.lower(dirX);
    auto iEnd   = localBox.upper(dirX);


    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_DOUBLE_EQ(srcFieldOld(ix), destField(ix));
    }
}
TEST_F(aFieldLinearTimeInterpolate, giveNewSrcForAlphaOne)
{
    double interpolateTime = 0.5;
    destNew->setTime(interpolateTime);

    auto& layout = srcNew->gridLayout;

    auto& srcFieldNew = srcNew->field;

    auto& destField = destNew->field;


    auto zero_transformation{zeroTransformation()};
    SAMRAI::hier::BoxContainer ghost_cntnr;
    FieldOverlap<dim> overlap{ghost_cntnr, zero_transformation};

    timeOp.timeInterpolate(*destNew, domain, overlap, *srcOld, *srcNew);


    bool const withGhost{true};
    auto box = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(domain, qty, layout,
                                                                          !withGhost);

    auto ghostBox = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(domain, qty, layout,
                                                                               withGhost);

    auto localBox = AMRToLocal(static_cast<std::add_const_t<decltype(box)>>(box), ghostBox);

    auto iStart = localBox.lower(dirX);
    auto iEnd   = localBox.upper(dirX);


    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_DOUBLE_EQ(srcFieldNew(ix), destField(ix));
    }
}

TEST_F(aFieldLinearTimeInterpolate, giveEvaluationOnTheInterpolateTimeForLinear)
{
    double interpolateTime = 0.2;
    destNew->setTime(interpolateTime);

    auto& layout = srcNew->gridLayout;


    auto& destField = destNew->field;

    auto zero_transformation{zeroTransformation()};
    SAMRAI::hier::BoxContainer ghost_cntnr;
    FieldOverlap<dim> overlap{ghost_cntnr, zero_transformation};

    timeOp.timeInterpolate(*destNew, domain, overlap, *srcOld, *srcNew);


    bool const withGhost{true};
    auto box = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(domain, qty, layout,
                                                                          !withGhost);

    auto ghostBox = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(domain, qty, layout,
                                                                               withGhost);

    auto localBox = AMRToLocal(static_cast<std::add_const_t<decltype(box)>>(box), ghostBox);

    auto iStart = localBox.lower(dirX);
    auto iEnd   = localBox.upper(dirX);


    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        auto position = layout.fieldNodeCoordinates(destField, origin, ix);

        EXPECT_DOUBLE_EQ(srcFunc(position[dirX], interpolateTime), destField(ix));
    }
}

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
