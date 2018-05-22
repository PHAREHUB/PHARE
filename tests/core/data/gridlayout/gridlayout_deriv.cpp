#include "gridlayout_deriv.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"


#include "gridlayout_deriv.h"

TEST_F(a1DDerivative, DXBY1D)
{
    auto layout        = GridLayout<GridLayoutImplYee<1, 1>>{{{0.1}}, {{50}}, Point<double, 1>{0.}};
    auto expDerValue   = read("dxBy_interpOrder_1.txt");
    auto psi_d         = layout.physicalStartIndex(QtyCentering::dual, Direction::X);
    auto pei_d         = layout.physicalEndIndex(QtyCentering::dual, Direction::X);
    constexpr auto dir = Direction::X;

    auto centering = layout.centering(By.physicalQuantity());
    if (centering[0] == QtyCentering::dual)
    {
        for (auto ix = psi_d; ix <= pei_d; ++ix)
        {
            auto localDerivative = layout.deriv(By, make_index(ix), DirectionTag<Direction::X>{});
            EXPECT_DOUBLE_EQ(expDerValue[ix], localDerivative);
        }
    }
}
