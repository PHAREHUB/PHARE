
#include <math.h>

#include "gridlayout_deriv.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_test.hpp"


std::vector<double> read(std::string filename)
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
    std::string filename = std::string("dxBy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");
    auto expDerValue            = read(filename);
    [[maybe_unused]] auto gei_d = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    for (auto ix = 0u; ix <= gei_d_X; ++ix)
    {
        Point<double, 1> point = this->layout.fieldNodeCoordinates(this->By, {0.}, ix);
        this->By(ix)           = std::cos(2 * M_PI / 5. * point[0]);
    }

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        auto localDerivative = this->layout.deriv(this->By, {ix}, DirectionTag<Direction::X>{});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[ix], 1e-12));
    }
}



TYPED_TEST(a1DDerivative, DXEZ1D)
{
    std::string filename = std::string("dxEz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);

    for (auto ix = 0u; ix <= gei_p_X; ++ix)
    {
        Point<double, 1> point = this->layout.fieldNodeCoordinates(this->Ez, {0.}, ix);
        this->Ez(ix)           = std::cos(2 * M_PI / 5. * point[0]);
    }

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        auto localDerivative = this->layout.deriv(this->Ez, {ix}, DirectionTag<Direction::X>{});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[ix], 1e-12));
    }
}



// -----------------------------------------------------------------------------
//              2D case
// -----------------------------------------------------------------------------

using layoutImpls2D
    = ::testing::Types<GridLayoutImplYee<2, 1>, GridLayoutImplYee<2, 2>, GridLayoutImplYee<2, 3>>;

TYPED_TEST_SUITE(a2DDerivative, layoutImpls2D);



TYPED_TEST(a2DDerivative, DXBY2D)
{
    std::string filename = std::string("dxBy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto expDerValue = read(filename);

    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

    for (auto ix = 0u; ix <= gei_d_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->By, {0., 0.}, ix, iy);
            this->By(ix, iy)
                = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    std::array<std::uint32_t, 2> nPts_
        = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::X);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->By, {ix, iy}, DirectionTag<Direction::X>{});
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
        }
    }
}


TYPED_TEST(a2DDerivative, DYBY2D)
{
    std::string filename = std::string("dyBy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto expDerValue = read(filename);

    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    for (auto ix = 0u; ix <= gei_d_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->By, {0., 0.}, ix, iy);
            this->By(ix, iy)
                = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Y);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->By, {ix, iy}, DirectionTag<Direction::Y>{});
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
        }
    }
}


TYPED_TEST(a2DDerivative, DXEZ2D)
{
    std::string filename = std::string("dxEz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);

    for (auto ix = 0u; ix <= gei_p_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            Point<double, 2> point = this->layout.fieldNodeCoordinates(this->Ez, {0., 0.}, ix, iy);
            this->Ez(ix, iy)
                = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::X);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->Ez, {ix, iy}, DirectionTag<Direction::X>{});
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
        }
    }
}


TYPED_TEST(a2DDerivative, DYEZ2D)
{
    std::string filename = std::string("dyEz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);

    for (auto ix = 0u; ix <= gei_p_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            auto point = this->layout.fieldNodeCoordinates(this->Ez, {0., 0.}, ix, iy);
            this->Ez(ix, iy)
                = std::cos(2 * M_PI / 5. * point[0]) * std::sin(2 * M_PI / 6. * point[1]);
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Y);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            auto localDerivative
                = this->layout.deriv(this->Ez, {ix, iy}, DirectionTag<Direction::Y>{});
            auto index_ = ix * nPts_[1] + iy;
            EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
        }
    }
}



// -----------------------------------------------------------------------------
//              3D case
// -----------------------------------------------------------------------------

using layoutImpls3D
    = ::testing::Types<GridLayoutImplYee<3, 1>, GridLayoutImplYee<3, 2>, GridLayoutImplYee<3, 3>>;

TYPED_TEST_SUITE(a3DDerivative, layoutImpls3D);



TYPED_TEST(a3DDerivative, DXBY3D)
{
    std::string filename = std::string("dxBy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto expDerValue = read(filename);

    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);
    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    auto pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (auto ix = 0u; ix <= gei_d_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = 0u; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(this->By, {0., 0., 0.}, ix, iy, iz);
                this->By(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                       * std::cos(2 * M_PI / 6. * point[1])
                                       * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::X);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->By, {ix, iy, iz}, DirectionTag<Direction::X>{});
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DDerivative, DYBY3D)
{
    std::string filename = std::string("dyBy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto expDerValue = read(filename);

    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);
    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);
    auto psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    auto pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (auto ix = 0u; ix <= gei_d_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = 0u; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(this->By, {0., 0., 0.}, ix, iy, iz);
                this->By(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                       * std::cos(2 * M_PI / 6. * point[1])
                                       * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Y);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->By, {ix, iy, iz}, DirectionTag<Direction::Y>{});
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DDerivative, DZBY3D)
{
    std::string filename = std::string("dzBy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto expDerValue = read(filename);

    auto gei_d_X = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);

    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_p_Z = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
    auto pei_p_Z = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Z);

    for (auto ix = 0u; ix <= gei_d_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = 0u; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(this->By, {0., 0., 0.}, ix, iy, iz);
                this->By(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                       * std::cos(2 * M_PI / 6. * point[1])
                                       * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::By, Direction::Z);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (auto iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->By, {ix, iy, iz}, DirectionTag<Direction::Z>{});
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DDerivative, DXEZ3D)
{
    std::string filename = std::string("dxEz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);
    auto psi_d_X = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d_X = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    auto pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (auto ix = 0u; ix <= gei_p_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = 0u; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Ez, {0., 0., 0.}, ix, iy, iz);
                this->Ez(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                       * std::cos(2 * M_PI / 6. * point[1])
                                       * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::X);

    for (auto ix = psi_d_X; ix <= pei_d_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->Ez, {ix, iy, iz}, DirectionTag<Direction::X>{});
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DDerivative, DYEZ3D)
{
    std::string filename = std::string("dyEz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);
    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_d_Y = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Y);
    auto pei_d_Y = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Y);
    auto psi_d_Z = this->layout.physicalStartIndex(QtyCentering::dual, Direction::Z);
    auto pei_d_Z = this->layout.physicalEndIndex(QtyCentering::dual, Direction::Z);

    for (auto ix = 0u; ix <= gei_p_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = 0u; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Ez, {0., 0., 0.}, ix, iy, iz);
                this->Ez(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                       * std::cos(2 * M_PI / 6. * point[1])
                                       * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Y);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_d_Y; iy <= pei_d_Y; ++iy)
        {
            for (auto iz = psi_d_Z; iz <= pei_d_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->Ez, {ix, iy, iz}, DirectionTag<Direction::Y>{});
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DDerivative, DZEZ3D)
{
    std::string filename = std::string("dzEz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");
    auto expDerValue = read(filename);

    auto gei_p_X = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto gei_p_Y = this->layout.ghostEndIndex(QtyCentering::primal, Direction::Y);
    auto gei_d_Z = this->layout.ghostEndIndex(QtyCentering::dual, Direction::Z);
    auto psi_p_X = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p_X = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);
    auto psi_p_Y = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Y);
    auto pei_p_Y = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Y);
    auto psi_p_Z = this->layout.physicalStartIndex(QtyCentering::primal, Direction::Z);
    auto pei_p_Z = this->layout.physicalEndIndex(QtyCentering::primal, Direction::Z);

    for (auto ix = 0u; ix <= gei_p_X; ++ix)
    {
        for (auto iy = 0u; iy <= gei_p_Y; ++iy)
        {
            for (auto iz = 0u; iz <= gei_d_Z; ++iz)
            {
                auto point = this->layout.fieldNodeCoordinates(this->Ez, {0., 0., 0.}, ix, iy, iz);
                this->Ez(ix, iy, iz) = std::sin(2 * M_PI / 5. * point[0])
                                       * std::cos(2 * M_PI / 6. * point[1])
                                       * std::sin(2 * M_PI / 12. * point[2]);
            }
        }
    }

    auto nPts_ = this->layout.allocSizeDerived(HybridQuantity::Scalar::Ez, Direction::Z);

    for (auto ix = psi_p_X; ix <= pei_p_X; ++ix)
    {
        for (auto iy = psi_p_Y; iy <= pei_p_Y; ++iy)
        {
            for (auto iz = psi_p_Z; iz <= pei_p_Z; ++iz)
            {
                auto localDerivative
                    = this->layout.deriv(this->Ez, {ix, iy, iz}, DirectionTag<Direction::Z>{});
                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;
                EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[index_], 1e-12));
            }
        }
    }
}
