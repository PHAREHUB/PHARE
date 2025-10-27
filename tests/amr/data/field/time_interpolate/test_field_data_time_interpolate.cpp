#include <type_traits>

#include "core/def/phare_mpi.hpp"

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "amr/resources_manager/amr_utils.hpp"



using namespace PHARE::core;
using namespace PHARE::amr;



template<typename T>
T srcFunc(T v)
{
    return v;
}

template<typename T, typename... Args>
T srcFunc(T v, Args... args)
{
    // constexpr std::size_t n = sizeof...(Args);
    // ASSERT_GT(n, 1);
    // ASSERT_LE(n, 3);
    return v + 10 * srcFunc(args...); // t + 10*x + 100*y + 1000*z
}


template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct aFieldLinearTimeInterpolate : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();

    using GridYee    = GridLayout<GridLayoutImplYee<dim, interp>>;
    using GridND     = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using FieldDataT = FieldData<GridYee, GridND>;

    FieldLinearTimeInterpolate<GridYee, GridND> timeOp{};
    HybridQuantity::Scalar qty{HybridQuantity::Scalar::Bx};
    SAMRAI::tbox::Dimension dimension{dim};
    static int countLocal;

    static constexpr auto lower = 0;
    static constexpr auto upper = 20;

    SAMRAI::hier::Box domain{SAMRAI::hier::Index{dimension, lower},
                             SAMRAI::hier::Index{dimension, upper}, SAMRAI::hier::BlockId{0},
                             SAMRAI::hier::LocalId{++countLocal},
                             SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank()};

    std::string const fieldName{"Bx"};
    SAMRAI::hier::IntVector ghost{dimension, GridYee::nbrGhosts()};
    static auto constexpr dl = ConstArray<double, dim>(0.01);

    static constexpr auto nbrCells = ConstArray<std::uint32_t, dim>(upper - lower + 1);

    Point<double, dim> origin{{0.}};

    std::shared_ptr<FieldDataT> srcOld;
    std::shared_ptr<FieldDataT> srcNew;
    std::shared_ptr<FieldDataT> destNew;


    aFieldLinearTimeInterpolate()
        : srcOld{std::make_shared<FieldDataT>(domain, ghost, fieldName, dl, nbrCells, origin, qty)}
        , srcNew{std::make_shared<FieldDataT>(domain, ghost, fieldName, dl, nbrCells, origin, qty)}
        , destNew{std::make_shared<FieldDataT>(domain, ghost, fieldName, dl, nbrCells, origin, qty)}
    {
        double const oldTime = 0.;
        double const newTime = 0.5;

        srcOld->setTime(oldTime);
        srcNew->setTime(newTime);

        auto& layout      = srcOld->gridLayout;
        auto& srcFieldOld = srcOld->field;
        auto& srcFieldNew = srcNew->field;

        for (auto const amr_idx : layout.AMRGhostBoxFor(qty))
        {
            auto const position = layout.fieldNodeCoordinates(srcFieldOld, amr_idx);
            auto const lcl_idx  = layout.AMRToLocal(amr_idx);

            if constexpr (dim == 1)
            {
                srcFieldOld(lcl_idx) = srcFunc(oldTime, position[dirX]);
                srcFieldNew(lcl_idx) = srcFunc(newTime, position[dirX]);
            }
            if constexpr (dim == 2)
            {
                srcFieldOld(lcl_idx) = srcFunc(oldTime, position[dirX], position[dirY]);
                srcFieldNew(lcl_idx) = srcFunc(newTime, position[dirX], position[dirY]);
            }
            if constexpr (dim == 3)
            {
                srcFieldOld(lcl_idx)
                    = srcFunc(oldTime, position[dirX], position[dirY], position[dirZ]);
                srcFieldNew(lcl_idx)
                    = srcFunc(newTime, position[dirX], position[dirY], position[dirZ]);
            }
        }
    }

    auto zeroTransformation() const
    {
        return SAMRAI::hier::Transformation(SAMRAI::hier::Transformation::NO_ROTATE,
                                            SAMRAI::hier::IntVector::getZero(dimension),
                                            SAMRAI::hier::BlockId(0), SAMRAI::hier::BlockId(0));
    }
};


using aFieldLinearTimeInterpolateInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<1>, InterpConst<2>>,
                     std::pair<DimConst<1>, InterpConst<3>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<2>, InterpConst<2>>, std::pair<DimConst<2>, InterpConst<3>>,
                     std::pair<DimConst<3>, InterpConst<1>>, std::pair<DimConst<3>, InterpConst<2>>,
                     std::pair<DimConst<3>, InterpConst<3>>>;

TYPED_TEST_SUITE(aFieldLinearTimeInterpolate, aFieldLinearTimeInterpolateInfos);


template<typename TypeInfo>
int aFieldLinearTimeInterpolate<TypeInfo>::countLocal = 0;



TYPED_TEST(aFieldLinearTimeInterpolate, giveOldSrcForAlphaZero)
{
    double const interpolateTime = 0.;
    this->destNew->setTime(interpolateTime);

    auto& layout      = this->srcOld->gridLayout;
    auto& srcFieldOld = this->srcOld->field;
    auto& destField   = this->destNew->field;

    auto zero_transformation{this->zeroTransformation()};
    SAMRAI::hier::BoxContainer ghost_cntnr;
    FieldOverlap overlap{ghost_cntnr, zero_transformation};

    this->timeOp.timeInterpolate(*(this->destNew), this->domain, overlap, *(this->srcOld),
                                 *(this->srcNew));


    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    auto box = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(this->domain, this->qty,
                                                                          layout);

    auto ghostBox_{this->domain};
    ghostBox_.grow(SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dim},
                                           static_cast<int>(GridYee::nbrGhosts())});
    auto ghostBox
        = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(ghostBox_, this->qty, layout);

    auto localBox = AMRToLocal(static_cast<std::add_const_t<decltype(box)>>(box), ghostBox);


    if constexpr (dim == 1)
    {
        auto iStartX = localBox.lower(dirX);
        auto iEndX   = localBox.upper(dirX);

        for (auto ix = iStartX; ix <= iEndX; ++ix)
        {
            EXPECT_DOUBLE_EQ(srcFieldOld(ix), destField(ix));
        }
    }
    if constexpr (dim == 2)
    {
        auto iStartX = localBox.lower(dirX);
        auto iEndX   = localBox.upper(dirX);
        auto iStartY = localBox.lower(dirY);
        auto iEndY   = localBox.upper(dirY);

        for (auto ix = iStartX; ix <= iEndX; ++ix)
        {
            for (auto iy = iStartY; iy <= iEndY; ++iy)
            {
                EXPECT_DOUBLE_EQ(srcFieldOld(ix, iy), destField(ix, iy));
            }
        }
    }
    if constexpr (dim == 3)
    {
        auto iStartX = localBox.lower(dirX);
        auto iEndX   = localBox.upper(dirX);
        auto iStartY = localBox.lower(dirY);
        auto iEndY   = localBox.upper(dirY);
        auto iStartZ = localBox.lower(dirZ);
        auto iEndZ   = localBox.upper(dirZ);

        for (auto ix = iStartX; ix <= iEndX; ++ix)
        {
            for (auto iy = iStartY; iy <= iEndY; ++iy)
            {
                for (auto iz = iStartZ; iz <= iEndZ; ++iz)
                {
                    EXPECT_DOUBLE_EQ(srcFieldOld(ix, iy, iz), destField(ix, iy, iz));
                }
            }
        }
    }
}



TYPED_TEST(aFieldLinearTimeInterpolate, giveNewSrcForAlphaOne)
{
    double interpolateTime = 0.5;
    this->destNew->setTime(interpolateTime);

    auto& layout      = this->srcOld->gridLayout;
    auto& srcFieldNew = this->srcNew->field;
    auto& destField   = this->destNew->field;

    auto zero_transformation{this->zeroTransformation()};
    SAMRAI::hier::BoxContainer ghost_cntnr;
    FieldOverlap overlap{ghost_cntnr, zero_transformation};

    this->timeOp.timeInterpolate(*(this->destNew), this->domain, overlap, *(this->srcOld),
                                 *(this->srcNew));


    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;

    auto box = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(this->domain, this->qty,
                                                                          layout);

    auto ghostBox_{this->domain};
    ghostBox_.grow(SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dim},
                                           static_cast<int>(GridYee::nbrGhosts())});
    auto ghostBox
        = FieldGeometry<GridYee, HybridQuantity::Scalar>::toFieldBox(ghostBox_, this->qty, layout);

    auto localBox = AMRToLocal(static_cast<std::add_const_t<decltype(box)>>(box), ghostBox);


    if constexpr (dim == 1)
    {
        auto iStartX = localBox.lower(dirX);
        auto iEndX   = localBox.upper(dirX);

        for (auto ix = iStartX; ix <= iEndX; ++ix)
        {
            EXPECT_DOUBLE_EQ(srcFieldNew(ix), destField(ix));
        }
    }
    if constexpr (dim == 2)
    {
        auto iStartX = localBox.lower(dirX);
        auto iEndX   = localBox.upper(dirX);
        auto iStartY = localBox.lower(dirY);
        auto iEndY   = localBox.upper(dirY);

        for (auto ix = iStartX; ix <= iEndX; ++ix)
        {
            for (auto iy = iStartY; iy <= iEndY; ++iy)
            {
                EXPECT_DOUBLE_EQ(srcFieldNew(ix, iy), destField(ix, iy));
            }
        }
    }
    if constexpr (dim == 3)
    {
        auto iStartX = localBox.lower(dirX);
        auto iEndX   = localBox.upper(dirX);
        auto iStartY = localBox.lower(dirY);
        auto iEndY   = localBox.upper(dirY);
        auto iStartZ = localBox.lower(dirZ);
        auto iEndZ   = localBox.upper(dirZ);

        for (auto ix = iStartX; ix <= iEndX; ++ix)
        {
            for (auto iy = iStartY; iy <= iEndY; ++iy)
            {
                for (auto iz = iStartZ; iz <= iEndZ; ++iz)
                {
                    EXPECT_DOUBLE_EQ(srcFieldNew(ix, iy, iz), destField(ix, iy, iz));
                }
            }
        }
    }
}


TYPED_TEST(aFieldLinearTimeInterpolate, giveEvaluationOnTheInterpolateTimeForLinear)
{
    static constexpr auto dim = typename TypeParam::first_type{}();

    double const interpolateTime = 0.2;
    this->destNew->setTime(interpolateTime);

    auto& layout    = this->srcOld->gridLayout;
    auto& destField = this->destNew->field;

    SAMRAI::hier::BoxContainer ghost_cntnr;
    FieldOverlap overlap{ghost_cntnr, this->zeroTransformation()};

    this->timeOp.timeInterpolate(*(this->destNew), this->domain, overlap, *(this->srcOld),
                                 *(this->srcNew));

    for (auto const [amr_idx, lcl_idx] : layout.domain_amr_lcl_idx(destField))
    {
        auto const position = layout.fieldNodeCoordinates(destField, amr_idx);

        if constexpr (dim == 1)
        {
            EXPECT_DOUBLE_EQ(srcFunc(interpolateTime, position[dirX]), destField(lcl_idx));
        }
        if constexpr (dim == 2)
        {
            EXPECT_DOUBLE_EQ(srcFunc(interpolateTime, position[dirX], position[dirY]),
                             destField(lcl_idx));
        }
        if constexpr (dim == 3)
        {
            EXPECT_DOUBLE_EQ(
                srcFunc(interpolateTime, position[dirX], position[dirY], position[dirZ]),
                destField(lcl_idx));
        }
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
