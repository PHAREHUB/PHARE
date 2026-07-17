
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
            for (std::size_t j = 0; j < rho.size(); ++j)
                x += rho.data()[j];
    }

    PHARE_LOG_LINE_SS("x " << x);

    // if (x != C * C * C * N)
    //     throw std::runtime_error("no");
}
