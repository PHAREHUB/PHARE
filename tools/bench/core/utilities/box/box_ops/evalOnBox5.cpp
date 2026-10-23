
#include "evalOnBox.hpp"


#include "mkn/avx/grid.hpp"
#include "mkn/avx/vector.hpp"


int main()
{
    PHARE_LOG_LINE_SS(__FILE__);


    using Array_t                      = mkn::avx::Vector_t<T>;
    auto constexpr static field_ghosts = GridLayout_t::nbrGhosts();
    static_assert(field_ghosts == 2 && "unexpected"); // unexpected otherwise

    GridLayout_t layout{C - 1};
    auto const dims = layout.allocSize(PHARE::core::HybridQuantity::Scalar::rho);
    auto const size = PHARE::core::product(dims);
    auto const sizes
        = PHARE::core::generate_from([&](auto const e) { return std::size_t{e}; }, dims);

    Array_t const rho(size, 1);
    Array_t tmp(size, 0);

    mkn::avx::Grid<T const, 3> const grid0{rho.data(), sizes};
    mkn::avx::Grid<T, 3> grid1{tmp.data(), sizes};
    auto const subgrid0 = grid0 >> field_ghosts;
    auto subgrid1       = grid1 >> field_ghosts;

    {
        PHARE::core::Timer t{N};
        for (auto i = 0u; i < N; ++i)
            subgrid1 += subgrid0;
    }

    auto const x = PHARE::core::sum(tmp);
    PHARE::core::abort_if_not(x == C * C * C * N && "fail");
}
