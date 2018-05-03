#include "gridlayout_deriv.h"
#include "gridlayout_params.h"
#include "gridlayout_test.h"

namespace PHARE
{
using GridLayoutTestYee1D = GridLayoutTest<Layout::Yee, 1, GridLayoutDerivParam>;

using GridLayoutDeriv1D = GridLayoutTestYee1D;


TEST_P(GridLayoutDeriv1D, DerivIsOK)
{
    for (std::size_t ix = 0; ix < param.expectedDeriv.size(); ++ix)
    {
        auto const& index = param.iCellDeriv[ix];
        EXPECT_DOUBLE_EQ((*param.derivedField)(index[0]), param.expectedDeriv[ix]);
    }
}


INSTANTIATE_TEST_CASE_P(DerivTest, GridLayoutDeriv1D,
                        ::testing::ValuesIn(createDerivParam<Layout::Yee, 1>()));

} // namespace PHARE
