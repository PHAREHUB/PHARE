
#include <math.h>

#include "gridlayout_deriv.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"

// -----------------------------------------------------------------------------
//              1D case
// -----------------------------------------------------------------------------

using layoutImpls1D
    = ::testing::Types<GridLayoutImplYee<1, 1>, GridLayoutImplYee<1, 2>, GridLayoutImplYee<1, 3>>;

TYPED_TEST_CASE(a1DDerivative, layoutImpls1D);



TYPED_TEST(a1DDerivative, DXBY1D)
{
    std::string  filename = std::string("dxBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_1d.txt");
    auto expDerValue = read(filename);
    auto gei_d       = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    uint32 gei_d_X       = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = 0; ix <= gei_d_X; ++ix)
    {
         Point<double, 1> point   = this->layout.fieldNodeCoordinates(this->By, Point<double, 1>{0.}, ix);
        this->By(ix) = std::cos(2 * M_PI / 5. * point[0]);
    }

    for (uint32 ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        auto localDerivative
            = this->layout.deriv(this->By, make_index(ix), DirectionTag<Direction::X>{});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[ix], 1e-12));
    }
}



TYPED_TEST(a1DDerivative, DXEZ1D)
{
    std::string  filename = std::string("dxEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_1d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X       = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

    auto psi_d_X       = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X       = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

    for (uint32 ix = 0; ix <= gei_p_X; ++ix)
    {
        Point<double, 1>  point   = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 1>{0.}, ix);
        this->Ez(ix) = std::cos(2 * M_PI / 5. * point[0]);
    }

    for (uint32 ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        auto localDerivative
            = this->layout.deriv(this->Ez, make_index(ix), DirectionTag<Direction::X>{});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[ix], 1e-12));
    }
}



// -----------------------------------------------------------------------------
//              2D case
// -----------------------------------------------------------------------------

using layoutImpls2D
    = ::testing::Types<GridLayoutImplYee<2, 1>, GridLayoutImplYee<2, 2>, GridLayoutImplYee<2, 3>>;

TYPED_TEST_CASE(a2DDerivative, layoutImpls2D);



TYPED_TEST(a2DDerivative, DXBY2D)
{
    std::string filename = std::string("dxBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_2d.txt");
    auto expDerValue = read(filename);

    uint32 gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

    uint32 psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    uint32 pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    uint32 psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    uint32 pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

    for (uint32 ix = 0; ix <= gei_d_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->By, Point<double, 2>{0., 0.}, ix, iy);
            this->By(ix, iy) = std::cos(2 * M_PI / 5. * point[0])*std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    std::array<uint32, 2> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::X);

    for (uint32 ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (uint32 iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->By, make_index(ix, iy), DirectionTag<Direction::X>{});
            uint32 index_= ix*nPts_[1]+iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));

        }
    }
}


TYPED_TEST(a2DDerivative, DYBY2D)
{
    std::string filename = std::string("dyBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_2d.txt");
    auto expDerValue = read(filename);

    uint32 gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

    uint32 psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    uint32 pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    uint32 psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    uint32 pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    for (uint32 ix = 0; ix <= gei_d_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->By, Point<double, 2>{0., 0.}, ix, iy);
            this->By(ix, iy) = std::cos(2 * M_PI / 5. * point[0])*std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    std::array<uint32, 2> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Y);

    for (uint32 ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (uint32 iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->By, make_index(ix, iy), DirectionTag<Direction::Y>{});
            uint32 index_= ix*nPts_[1]+iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
        }
    }
}


TYPED_TEST(a2DDerivative, DXEZ2D)
{
    std::string filename = std::string("dxEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_2d.txt");
    auto expDerValue = read(filename);

    uint32 gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

    uint32 psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    uint32 pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    uint32 psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    uint32 pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

    for (uint32 ix = 0; ix <= gei_p_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 2>{0., 0.}, ix, iy);
            this->Ez(ix, iy) = std::cos(2 * M_PI / 5. * point[0])*std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    std::array<uint32, 2> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::X);

    for (uint32 ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (uint32 iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->Ez, make_index(ix, iy), DirectionTag<Direction::X>{});
            uint32 index_= ix*nPts_[1]+iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));

        }
    }
}


TYPED_TEST(a2DDerivative, DYEZ2D)
{
    std::string filename = std::string("dyEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_2d.txt");
    auto expDerValue = read(filename);

    uint32 gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

    uint32 psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    uint32 pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    uint32 psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    uint32 pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    for (uint32 ix = 0; ix <= gei_p_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 2>{0., 0.}, ix, iy);
            this->Ez(ix, iy) = std::cos(2 * M_PI / 5. * point[0])*std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    std::array<uint32, 2> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Y);

    for (uint32 ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (uint32 iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->Ez, make_index(ix, iy), DirectionTag<Direction::Y>{});
            uint32 index_= ix*nPts_[1]+iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));

        }
    }
}



// -----------------------------------------------------------------------------
//              3D case
// -----------------------------------------------------------------------------

using layoutImpls3D
    = ::testing::Types<GridLayoutImplYee<3, 1>, GridLayoutImplYee<3, 2>, GridLayoutImplYee<3, 3>>;

TYPED_TEST_CASE(a3DDerivative, layoutImpls3D);



TYPED_TEST(a3DDerivative, DXBY3D)
{
    std::string filename = std::string("dxBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_3d.txt");
    auto expDerValue = read(filename);

    uint32 gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    uint32 gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    uint32 psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    uint32 pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    uint32 psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    uint32 pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    uint32 psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    uint32 pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (uint32 ix = 0; ix <= gei_d_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            for (uint32 iz = 0; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(this->By, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                this->By(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])*std::cos(2 * M_PI / 6. * point[1])*std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    std::array<uint32, 3> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::X);

    for (uint32 ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (uint32 iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (uint32 iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->By, make_index(ix, iy, iz), DirectionTag<Direction::X>{});
                uint32 index_= ix*nPts_[1]*nPts_[2]+iy*nPts_[2]+iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }

        }
    }
}



TYPED_TEST(a3DDerivative, DYBY3D)
{
    std::string filename = std::string("dyBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_3d.txt");
    auto expDerValue = read(filename);

    uint32 gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    uint32 gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    uint32 psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    uint32 pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    uint32 psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    uint32 pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);
    uint32 psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    uint32 pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (uint32 ix = 0; ix <= gei_d_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            for (uint32 iz = 0; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(this->By, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                this->By(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])*std::cos(2 * M_PI / 6. * point[1])*std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    std::array<uint32, 3> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Y);

    for (uint32 ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (uint32 iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (uint32 iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->By, make_index(ix, iy, iz), DirectionTag<Direction::Y>{});
                uint32 index_= ix*nPts_[1]*nPts_[2]+iy*nPts_[2]+iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }

        }
    }
}



TYPED_TEST(a3DDerivative, DZBY3D)
{
    std::string filename = std::string("dzBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_3d.txt");
    auto expDerValue = read(filename);

    uint32 gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    uint32 gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    uint32 psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    uint32 pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    uint32 psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    uint32 pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    uint32 psi_p_Z = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
    uint32 pei_p_Z = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Z);

    for (uint32 ix = 0; ix <= gei_d_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            for (uint32 iz = 0; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(this->By, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                this->By(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])*std::cos(2 * M_PI / 6. * point[1])*std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    std::array<uint32, 3> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Z);

    for (uint32 ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (uint32 iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (uint32 iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->By, make_index(ix, iy, iz), DirectionTag<Direction::Z>{});
                uint32 index_= ix*nPts_[1]*nPts_[2]+iy*nPts_[2]+iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }

        }
    }
}



TYPED_TEST(a3DDerivative, DXEZ3D)
{
    std::string filename = std::string("dxEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_3d.txt");
    auto expDerValue = read(filename);

    uint32 gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    uint32 gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    uint32 psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    uint32 pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    uint32 psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    uint32 pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    uint32 psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    uint32 pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (uint32 ix = 0; ix <= gei_p_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            for (uint32 iz = 0; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                this->Ez(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])*std::cos(2 * M_PI / 6. * point[1])*std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    std::array<uint32, 3> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::X);

    for (uint32 ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (uint32 iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (uint32 iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->Ez, make_index(ix, iy, iz), DirectionTag<Direction::X>{});
                uint32 index_= ix*nPts_[1]*nPts_[2]+iy*nPts_[2]+iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }

        }
    }
}



TYPED_TEST(a3DDerivative, DYEZ3D)
{
    std::string filename = std::string("dyEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_3d.txt");
    auto expDerValue = read(filename);

    uint32 gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    uint32 gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    uint32 psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    uint32 pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    uint32 psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    uint32 pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);
    uint32 psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    uint32 pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (uint32 ix = 0; ix <= gei_p_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            for (uint32 iz = 0; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                this->Ez(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])*std::cos(2 * M_PI / 6. * point[1])*std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    std::array<uint32, 3> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Y);

    for (uint32 ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (uint32 iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (uint32 iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->Ez, make_index(ix, iy, iz), DirectionTag<Direction::Y>{});
                uint32 index_= ix*nPts_[1]*nPts_[2]+iy*nPts_[2]+iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }

        }
    }
}



TYPED_TEST(a3DDerivative, DZEZ3D)
{
    std::string filename = std::string("dzEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string("_3d.txt");
    auto expDerValue = read(filename);

    uint32 gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    uint32 gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    uint32 gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    uint32 psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    uint32 pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    uint32 psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    uint32 pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    uint32 psi_p_Z = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
    uint32 pei_p_Z = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Z);

    for (uint32 ix = 0; ix <= gei_p_X; ++ix)
    {
        for (uint32 iy = 0; iy <= gei_p_Y; ++iy)
        {
            for (uint32 iz = 0; iz <= gei_d_Z; ++iz)
            {
                Point<double, 3> point = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 3>{0., 0., 0.}, ix, iy, iz);
                this->Ez(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])*std::cos(2 * M_PI / 6. * point[1])*std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    std::array<uint32, 3> nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Z);

    for (uint32 ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (uint32 iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (uint32 iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->Ez, make_index(ix, iy, iz), DirectionTag<Direction::Z>{});
                uint32 index_= ix*nPts_[1]*nPts_[2]+iy*nPts_[2]+iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }

        }
    }
}

