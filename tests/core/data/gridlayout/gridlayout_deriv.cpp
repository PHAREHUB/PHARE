

#include "core/data/grid/gridlayoutdefs.hpp"

#include "gridlayout_deriv.hpp"

#include <math.h>
#include <string>


std::vector<double> read(std::string const filename)
{
    std::ifstream readFile(filename);
    assert(readFile.is_open());
    std::vector<double> x;

    std::copy(std::istream_iterator<double>(readFile), std::istream_iterator<double>(),
              std::back_inserter(x));
    return x;
}


// -----------------------------------------------------------------------------
//              1D case
// -----------------------------------------------------------------------------

using layoutImpls1D
    = ::testing::Types<GridLayoutImplYee<1, 1>, GridLayoutImplYee<1, 2>, GridLayoutImplYee<1, 3>>;

TYPED_TEST_SUITE(a1DDerivative, layoutImpls1D);



TYPED_TEST(a1DDerivative, DXBY1D)
{
    auto const filename = std::string("dxBy_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");
    auto const expDerValue = read(filename);


    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->By))
    {
        auto const point  = this->layout.fieldNodeCoordinates(this->By, amr_idx);
        this->By(lcl_idx) = std::cos(2 * M_PI / 5. * point[0]);
    }

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::X>(this->By, {lcl_idx});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[lcl_idx[0]], 1e-12));
    }
}



TYPED_TEST(a1DDerivative, DXEZ1D)
{
    auto const filename = std::string("dxEz_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");
    auto const expDerValue = read(filename);

    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->Ez))
    {
        auto const point  = this->layout.fieldNodeCoordinates(this->Ez, amr_idx);
        this->Ez(lcl_idx) = std::cos(2 * M_PI / 5. * point[0]);
    }

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::X>(this->Ez, {lcl_idx});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[lcl_idx[0]], 1e-12));
    }
}



// // -----------------------------------------------------------------------------
// //              2D case
// // -----------------------------------------------------------------------------

using layoutImpls2D
    = ::testing::Types<GridLayoutImplYee<2, 1>, GridLayoutImplYee<2, 2>, GridLayoutImplYee<2, 3>>;

TYPED_TEST_SUITE(a2DDerivative, layoutImpls2D);



TYPED_TEST(a2DDerivative, DXBY2D)
{
    auto const filename = std::string("dxBy_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto const expDerValue = read(filename);


    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->By))
    {
        Point<double, 2> point = this->layout.fieldNodeCoordinates(this->By, amr_idx);
        this->By(lcl_idx) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::X);


    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::X>(this->By, lcl_idx);
        auto index_                = lcl_idx[0] * nPts_[1] + lcl_idx[1];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}


TYPED_TEST(a2DDerivative, DYBY2D)
{
    auto const filename = std::string("dyBy_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto const expDerValue = read(filename);


    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->By))
    {
        Point<double, 2> point = this->layout.fieldNodeCoordinates(this->By, amr_idx);
        this->By(lcl_idx) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Y);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::Y>(this->By, lcl_idx);
        auto index_                = lcl_idx[0] * nPts_[1] + lcl_idx[1];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}


TYPED_TEST(a2DDerivative, DXEZ2D)
{
    auto const filename = std::string("dxEz_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto const expDerValue = read(filename);


    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->Ez))
    {
        Point<double, 2> point = this->layout.fieldNodeCoordinates(this->Ez, amr_idx);
        this->Ez(lcl_idx) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::X);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::X>(this->Ez, lcl_idx);
        auto index_                = lcl_idx[0] * nPts_[1] + lcl_idx[1];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}


TYPED_TEST(a2DDerivative, DYEZ2D)
{
    auto const filename = std::string("dyEz_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto const expDerValue = read(filename);


    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->Ez))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->Ez, amr_idx);
        this->Ez(lcl_idx) = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Y);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::Y>(this->Ez, lcl_idx);
        auto index_                = lcl_idx[0] * nPts_[1] + lcl_idx[1];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}



// // -----------------------------------------------------------------------------
// //              3D case
// // -----------------------------------------------------------------------------

using layoutImpls3D
    = ::testing::Types<GridLayoutImplYee<3, 1>, GridLayoutImplYee<3, 2>, GridLayoutImplYee<3, 3>>;

TYPED_TEST_SUITE(a3DDerivative, layoutImpls3D);



TYPED_TEST(a3DDerivative, DXBY3D)
{
    auto const filename = std::string("dxBy_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto const expDerValue = read(filename);


    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->By))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->By, amr_idx);
        this->By(lcl_idx) = std::sin(2 * M_PI / 5. * point[0]) * std::cos(2 * M_PI / 6. * point[1])
                            * std::sin(2 * M_PI / 12. * point[2]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::X);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::X>(this->By, lcl_idx);
        auto index_ = lcl_idx[0] * nPts_[1] * nPts_[2] + lcl_idx[1] * nPts_[2] + lcl_idx[2];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}



TYPED_TEST(a3DDerivative, DYBY3D)
{
    auto const filename = std::string("dyBy_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto const expDerValue = read(filename);

    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->By))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->By, amr_idx);
        this->By(lcl_idx) = std::sin(2 * M_PI / 5. * point[0]) * std::cos(2 * M_PI / 6. * point[1])
                            * std::sin(2 * M_PI / 12. * point[2]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Y);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::Y>(this->By, lcl_idx);
        auto index_ = lcl_idx[0] * nPts_[1] * nPts_[2] + lcl_idx[1] * nPts_[2] + lcl_idx[2];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}



TYPED_TEST(a3DDerivative, DZBY3D)
{
    auto const filename = std::string("dzBy_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto const expDerValue = read(filename);

    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->By))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->By, amr_idx);
        this->By(lcl_idx) = std::sin(2 * M_PI / 5. * point[0]) * std::cos(2 * M_PI / 6. * point[1])
                            * std::sin(2 * M_PI / 12. * point[2]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Z);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::Z>(this->By, lcl_idx);
        auto index_ = lcl_idx[0] * nPts_[1] * nPts_[2] + lcl_idx[1] * nPts_[2] + lcl_idx[2];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}



TYPED_TEST(a3DDerivative, DXEZ3D)
{
    auto const filename = std::string("dxEz_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto const expDerValue = read(filename);

    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->Ez))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->Ez, amr_idx);
        this->Ez(lcl_idx) = std::sin(2 * M_PI / 5. * point[0]) * std::cos(2 * M_PI / 6. * point[1])
                            * std::sin(2 * M_PI / 12. * point[2]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::X);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::X>(this->Ez, lcl_idx);
        auto index_ = lcl_idx[0] * nPts_[1] * nPts_[2] + lcl_idx[1] * nPts_[2] + lcl_idx[2];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}



TYPED_TEST(a3DDerivative, DYEZ3D)
{
    auto const filename = std::string("dyEz_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto const expDerValue = read(filename);

    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->Ez))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->Ez, amr_idx);
        this->Ez(lcl_idx) = std::sin(2 * M_PI / 5. * point[0]) * std::cos(2 * M_PI / 6. * point[1])
                            * std::sin(2 * M_PI / 12. * point[2]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Y);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::Y>(this->Ez, lcl_idx);
        auto index_ = lcl_idx[0] * nPts_[1] * nPts_[2] + lcl_idx[1] * nPts_[2] + lcl_idx[2];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}



TYPED_TEST(a3DDerivative, DZEZ3D)
{
    auto const filename = std::string("dzEz_interpOrder_")
                          + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto const expDerValue = read(filename);

    for (auto const [amr_idx, lcl_idx] : this->layout.ghost_amr_lcl_idx(this->Ez))
    {
        auto point        = this->layout.fieldNodeCoordinates(this->Ez, amr_idx);
        this->Ez(lcl_idx) = std::sin(2 * M_PI / 5. * point[0]) * std::cos(2 * M_PI / 6. * point[1])
                            * std::sin(2 * M_PI / 12. * point[2]);
    }

    auto const nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Z);

    for (auto const lcl_idx : this->layout.AMRToLocal(this->layout.AMRBox()))
    {
        auto const localDerivative = this->layout.template deriv<Direction::Z>(this->Ez, lcl_idx);
        auto index_ = lcl_idx[0] * nPts_[1] * nPts_[2] + lcl_idx[1] * nPts_[2] + lcl_idx[2];
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
    }
}
