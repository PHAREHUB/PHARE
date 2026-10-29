
#include "evalOnBox.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    using namespace PHARE::core;

    GridLayout_t layout{C - 1};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    Grid_t tmp{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 0};

    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
        {
            auto const rslabs = BoxSpan{layout.domainBoxFor(rho), rho};
            auto tslabs       = BoxSpan{layout.domainBoxFor(rho), tmp};
            auto tslabit      = tslabs.begin();

            for (auto rslabit = rslabs.begin(); rslabit != rslabs.end(); ++rslabit, ++tslabit)
            {
                auto const rrow = *rslabit;
                auto trow       = *tslabit;
                auto trowit     = trow.begin();

                for (auto rrowit = rrow.begin(); rrowit != rrow.end(); ++rrowit, ++trowit)
                    for (std::size_t i = 0; i < (*rrowit).size(); ++i)
                        (*trowit)[i] += (*rrowit)[i];
            }
        }
    }

    auto const x = PHARE::core::sum(tmp);
    PHARE::core::abort_if_not(x == C * C * C * N && "fail");
}
