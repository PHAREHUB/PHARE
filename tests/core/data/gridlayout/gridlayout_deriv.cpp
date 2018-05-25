
#include <math.h>

#include "gridlayout_deriv.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"

using layoutImpls
    = ::testing::Types<GridLayoutImplYee<1, 1>, GridLayoutImplYee<1, 2>, GridLayoutImplYee<1, 3>>;

TYPED_TEST_CASE(a1DDerivative, layoutImpls);



TYPED_TEST(a1DDerivative, DXBY1D)
{
    auto filename = std::string("dxBy_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string(".txt");
    auto expDerValue = read(filename);
    auto psi_d       = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d       = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto gei_d       = this->layout.ghostEndIndex(QtyCentering::dual, Direction::X);

    auto psi_p = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    constexpr auto dir = Direction::X;

    for (auto ix = 0; ix <= gei_d; ++ix)
    {
        auto point   = this->layout.fieldNodeCoordinates(this->By, Point<double, 1>{0.}, ix);
        this->By(ix) = std::cos(2 * M_PI / 5. * point[0]);
    }


    for (auto ix = psi_p; ix <= pei_p; ++ix)
    {
        auto localDerivative
            = this->layout.deriv(this->By, make_index(ix), DirectionTag<Direction::X>{});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[ix], 1e-12));
    }
}



TYPED_TEST(a1DDerivative, DXEZ1D)
{
    auto filename = std::string("dxEz_interpOrder_") + std::to_string(TestFixture::interp_order)
                    + std::string(".txt");
    auto expDerValue = read(filename);
    auto psi_d       = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d       = this->layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    auto gei_p       = this->layout.ghostEndIndex(QtyCentering::primal, Direction::X);
    auto psi_p       = this->layout.physicalStartIndex(QtyCentering::primal, Direction::X);
    auto pei_p       = this->layout.physicalEndIndex(QtyCentering::primal, Direction::X);

    constexpr auto dir = Direction::X;

    for (auto ix = 0; ix <= gei_p; ++ix)
    {
        auto point   = this->layout.fieldNodeCoordinates(this->Ez, Point<double, 1>{0.}, ix);
        this->Ez(ix) = std::cos(2 * M_PI / 5. * point[0]);
    }


    for (auto ix = psi_d; ix <= pei_d; ++ix)
    {
        auto localDerivative
            = this->layout.deriv(this->Ez, make_index(ix), DirectionTag<Direction::X>{});
        EXPECT_THAT(localDerivative, ::testing::DoubleNear(expDerValue[ix], 1e-12));
    }
}
