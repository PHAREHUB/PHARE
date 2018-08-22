
// -------------------------------------
//     FieldLinearTimeInterpolate
// -------------------------------------

#include "data/field/field_data.h"
#include "data/field/field_geometry.h"

#include <SAMRAI/hier/TimeInterpolateOperator.h>


namespace PHARE
{
template<typename GridLayoutImpl, typename FieldT>
class FieldLinearTimeInterpolate : public SAMRAI::hier::TimeInterpolateOperator
{
public:
    FieldLinearTimeInterpolate()
        : SAMRAI::hier::TimeInterpolateOperator{"FieldLinearTimeInterpolate"}
    {
    }




    virtual ~FieldLinearTimeInterpolate() = default;




    virtual void timeInterpolate(SAMRAI::hier::PatchData& destData, SAMRAI::hier::Box const& where,
                                 SAMRAI::hier::PatchData const& srcDataOld,
                                 SAMRAI::hier::PatchData const& srcDataNew) const override
    {
        //
    }




private:
    static std::size_t constexpr dim = GridLayoutImpl::dimension;

    using PhysicalQuantity = decltype(std::declval<FieldT>().physicalQuantity());
    using FieldDataT       = FieldData<GridLayoutImpl, FieldT>;
};



} // namespace PHARE


// -------------------------------------
//   FieldLinearTimeInterpolate test
// -------------------------------------

#include <type_traits>

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


// #include "amr/data/field/time_interpolate/field_linear_time_interpolate.h"

#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "hybrid/hybrid_quantities.h"
#include "tools/amr_utils.h"



using namespace PHARE;

std::size_t constexpr dim         = 1;
std::size_t constexpr interpOrder = 1;


using GridYee = GridLayoutImplYee<dim, interpOrder>;

using Field1D = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

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

    std::array<uint32, dim> nbrCells{{static_cast<uint32>(domain.numberCells(dirX))}};

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
};

int aFieldLinearTimeInterpolate::countLocal = 0;

TEST_F(aFieldLinearTimeInterpolate, giveOldSrcForAlphaZero)
{
    double interpolateTime = 0.;
    destNew->setTime(interpolateTime);

    auto& layout = srcOld->gridLayout;


    auto& srcFieldOld = srcOld->field;

    auto& destField = destNew->field;


    timeOp.timeInterpolate(*destNew, domain, *srcOld, *srcNew);


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


    timeOp.timeInterpolate(*destNew, domain, *srcOld, *srcNew);


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
