
#include <math.h>

#include "gridlayout_laplacian.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_test.hpp"


// -----------------------------------------------------------------------------
//              1D case
// -----------------------------------------------------------------------------

using layoutImpls1D
    = ::testing::Types<GridLayoutImplYee<1, 1>, GridLayoutImplYee<1, 2>, GridLayoutImplYee<1, 3>>;

TYPED_TEST_SUITE(a1DLaplacian, layoutImpls1D);


TYPED_TEST(a1DLaplacian, LaplacianJx1D)
{
    std::string filename = std::string("lapJx_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jx, Direction::X);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jx, Direction::X);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        auto point   = this->layout.fieldNodeCoordinates(this->Jx, Point{0.}, ix);
        this->Jx(ix) = std::sinh(0.1 * point[0]);
    }

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        auto localLaplacian = this->layout.laplacian(this->Jx, make_index(ix));
        EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[ix], 1e-12));
    }
}



TYPED_TEST(a1DLaplacian, LaplacianJy1D)
{
    std::string filename = std::string("lapJy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::X);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::X);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        auto point   = this->layout.fieldNodeCoordinates(this->Jy, Point{0.}, ix);
        this->Jy(ix) = std::sinh(0.3 * point[0]);
    }

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        auto localLaplacian = this->layout.laplacian(this->Jy, make_index(ix));
        EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[ix], 1e-12));
    }
}



TYPED_TEST(a1DLaplacian, LaplacianJz1D)
{
    std::string filename = std::string("lapJz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_1d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jz, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jz, Direction::X);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jz, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jz, Direction::X);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        auto point   = this->layout.fieldNodeCoordinates(this->Jz, Point{0.}, ix);
        this->Jz(ix) = std::sinh(0.2 * point[0]);
    }

    for (std::uint32_t ix = psi_X; ix <= pei_X; ++ix)
    {
        auto localLaplacian = this->layout.laplacian(this->Jz, make_index(ix));
        EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[ix], 1e-12));
    }
}



// -----------------------------------------------------------------------------
//              2D case
// -----------------------------------------------------------------------------

using layoutImpls2D
    = ::testing::Types<GridLayoutImplYee<2, 1>, GridLayoutImplYee<2, 2>, GridLayoutImplYee<2, 3>>;

TYPED_TEST_SUITE(a2DLaplacian, layoutImpls2D);


TYPED_TEST(a2DLaplacian, LaplacianJx2D)
{
    std::string filename = std::string("lapJx_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto gsi_Y = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jx, Direction::Y);
    auto gei_Y = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jx, Direction::Y);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jx, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jx, Direction::Y);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
        {
            auto point       = this->layout.fieldNodeCoordinates(this->Jx, Point{0., 0.}, ix, iy);
            this->Jx(ix, iy) = std::sinh(0.1 * point[0]) * std::cosh(0.1 * point[1]);
        }
    }

    std::array<std::uint32_t, 2> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jx);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            auto localLaplacian = this->layout.laplacian(this->Jx, make_index(ix, iy));

            auto index_ = ix * nPts_[1] + iy;

            EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[index_], 1e-12));
        }
    }
}


TYPED_TEST(a2DLaplacian, LaplacianJy2D)
{
    std::string filename = std::string("lapJy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");

    auto expLapValue = read(filename);

    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gsi_Y = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto gei_Y = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::Y);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::Y);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
        {
            auto point       = this->layout.fieldNodeCoordinates(this->Jy, Point{0., 0.}, ix, iy);
            this->Jy(ix, iy) = std::sinh(0.3 * point[0]) * std::cosh(0.3 * point[1]);
        }
    }

    std::array<std::uint32_t, 2> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jy);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            auto localLaplacian = this->layout.laplacian(this->Jy, make_index(ix, iy));

            auto index_ = ix * nPts_[1] + iy;

            EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[index_], 1e-12));
        }
    }
}


TYPED_TEST(a2DLaplacian, LaplacianJz2D)
{
    std::string filename = std::string("lapJz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_2d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gsi_Y = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto gei_Y = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::Y);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::Y);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
        {
            auto point       = this->layout.fieldNodeCoordinates(this->Jz, Point{0., 0.}, ix, iy);
            this->Jz(ix, iy) = std::sinh(0.2 * point[0]) * std::cosh(0.2 * point[1]);
        }
    }

    std::array<std::uint32_t, 2> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jz);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            auto localLaplacian = this->layout.laplacian(this->Jz, make_index(ix, iy));

            auto index_ = ix * nPts_[1] + iy;

            EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[index_], 1e-12));
        }
    }
}




// -----------------------------------------------------------------------------
//              3D case
// -----------------------------------------------------------------------------

using layoutImpls3D
    = ::testing::Types<GridLayoutImplYee<3, 1>, GridLayoutImplYee<3, 2>, GridLayoutImplYee<3, 3>>;

TYPED_TEST_SUITE(a3DLaplacian, layoutImpls3D);


TYPED_TEST(a3DLaplacian, LaplacianJx3D)
{
    std::string filename = std::string("lapJx_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto gsi_Y = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jx, Direction::Y);
    auto gei_Y = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jx, Direction::Y);
    auto gsi_Z = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jx, Direction::Z);
    auto gei_Z = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jx, Direction::Z);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jx, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jx, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jx, Direction::Y);
    auto psi_Z = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jx, Direction::Z);
    auto pei_Z = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jx, Direction::Z);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
        {
            for (auto iz = gsi_Z; iz <= gei_Z; ++iz)
            {
                auto point
                    = this->layout.fieldNodeCoordinates(this->Jx, Point{0., 0., 0.}, ix, iy, iz);
                this->Jx(ix, iy, iz) = std::sinh(0.1 * point[0]) * std::cosh(0.1 * point[1])
                                       * std::tanh(0.1 * point[2]);
            }
        }
    }

    std::array<std::uint32_t, 3> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jx);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
            {
                auto localLaplacian = this->layout.laplacian(this->Jx, make_index(ix, iy, iz));

                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;

                EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DLaplacian, LaplacianJy3D)
{
    std::string filename = std::string("lapJy_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto gsi_Y = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto gei_Y = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto gsi_Z = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jy, Direction::Z);
    auto gei_Z = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jy, Direction::Z);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::Y);
    auto psi_Z = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jy, Direction::Z);
    auto pei_Z = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jy, Direction::Z);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
        {
            for (auto iz = gsi_Z; iz <= gei_Z; ++iz)
            {
                auto point
                    = this->layout.fieldNodeCoordinates(this->Jy, Point{0., 0., 0.}, ix, iy, iz);
                this->Jy(ix, iy, iz) = std::sinh(0.3 * point[0]) * std::cosh(0.3 * point[1])
                                       * std::tanh(0.3 * point[2]);
            }
        }
    }

    std::array<std::uint32_t, 3> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jy);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
            {
                auto localLaplacian = this->layout.laplacian(this->Jy, make_index(ix, iy, iz));

                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;

                EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[index_], 1e-12));
            }
        }
    }
}



TYPED_TEST(a3DLaplacian, LaplacianJz3D)
{
    std::string filename = std::string("lapJz_interpOrder_")
                           + std::to_string(TestFixture::interp_order) + std::string("_3d.txt");

    auto expLapValue = read(filename);

    auto gsi_X = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jz, Direction::X);
    auto gei_X = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jz, Direction::X);
    auto gsi_Y = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jz, Direction::Y);
    auto gei_Y = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jz, Direction::Y);
    auto gsi_Z = this->layout.ghostStartIndex(HybridQuantity::Scalar::Jz, Direction::Z);
    auto gei_Z = this->layout.ghostEndIndex(HybridQuantity::Scalar::Jz, Direction::Z);

    auto psi_X = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jz, Direction::X);
    auto pei_X = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jz, Direction::X);
    auto psi_Y = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jz, Direction::Y);
    auto pei_Y = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jz, Direction::Y);
    auto psi_Z = this->layout.physicalStartIndex(HybridQuantity::Scalar::Jz, Direction::Z);
    auto pei_Z = this->layout.physicalEndIndex(HybridQuantity::Scalar::Jz, Direction::Z);

    for (auto ix = gsi_X; ix <= gei_X; ++ix)
    {
        for (auto iy = gsi_Y; iy <= gei_Y; ++iy)
        {
            for (auto iz = gsi_Z; iz <= gei_Z; ++iz)
            {
                auto point
                    = this->layout.fieldNodeCoordinates(this->Jz, Point{0., 0., 0.}, ix, iy, iz);
                this->Jz(ix, iy, iz) = std::sinh(0.2 * point[0]) * std::cosh(0.2 * point[1])
                                       * std::tanh(0.2 * point[2]);
            }
        }
    }

    std::array<std::uint32_t, 3> nPts_ = this->layout.allocSize(HybridQuantity::Scalar::Jz);

    for (auto ix = psi_X; ix <= pei_X; ++ix)
    {
        for (auto iy = psi_Y; iy <= pei_Y; ++iy)
        {
            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
            {
                auto localLaplacian = this->layout.laplacian(this->Jz, make_index(ix, iy, iz));

                auto index_ = ix * nPts_[1] * nPts_[2] + iy * nPts_[2] + iz;

                EXPECT_THAT(localLaplacian, ::testing::DoubleNear(expLapValue[index_], 1e-12));
            }
        }
    }
}
