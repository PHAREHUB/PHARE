
#include "evalOnBox.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    using namespace PHARE::core;

    GridLayout_t layout{C - 1};
    auto const local_box
        = PHARE::core::box_from_zero_to_upper(*layout.AMRBox().shape().as_unsigned());
    std::array<std::uint32_t, dim> const dims = local_box.shape() + 1;
    Grid_t rho{"rho", PHARE::core::HybridQuantity::Scalar::rho, dims};
    Grid_t tmp{"rho", PHARE::core::HybridQuantity::Scalar::rho, dims};

    rho.fill(1);
    tmp.fill(0);

    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
        {
            auto const rslabs   = BoxSpan{local_box, rho};
            auto tslabs         = BoxSpan{local_box, tmp};
            auto tslabit        = tslabs.begin();
            auto const rslabend = rslabs.end();

            for (auto rslabit = rslabs.begin(); rslabit != rslabend; ++rslabit, ++tslabit)
            {
                auto const rrow    = *rslabit;
                auto trow          = *tslabit;
                auto trowit        = trow.begin();
                auto const rrowend = rrow.end();

                for (auto rrowit = rrow.begin(); rrowit != rrowend; ++rrowit, ++trowit)
                    for (std::size_t i = 0; i < (*rrowit).size(); ++i)
                        (*trowit)[i] += (*rrowit)[i];
            }
        }
    }

    auto const x = PHARE::core::sum(tmp);
    PHARE::core::abort_if_not(x == C * C * C * N && "fail");
}
