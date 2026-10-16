
#include "evalOnBox.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    GridLayout_t const layout{C - 1};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};

    T x = 0;
    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            layout.evalOnBox(rho, [&](auto const&... args) { x += rho(args...); });
    }

    if (x != C * C * C * N)
        throw std::runtime_error("no");
}
