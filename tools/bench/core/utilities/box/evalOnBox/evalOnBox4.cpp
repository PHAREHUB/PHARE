
#include "evalOnBox.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);

    using namespace PHARE::core;

    GridLayout_t layout{C - 1};
    Grid_t rho{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 1};
    T x = 0;
    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            for (auto const& slab : BoxSpan{layout.domainBoxFor(rho), rho})
                for (auto const& row : slab)
                    for (auto const& r : row)
                        x += r;
    }

    if (x != C * C * C * N)
        throw std::runtime_error("no");
}
