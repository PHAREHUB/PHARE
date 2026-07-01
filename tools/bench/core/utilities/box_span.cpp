#include "bench_box.hpp"



namespace PHARE::bench::core
{

int test()
{
    std::size_t constexpr static FOR = 10;

    GridLayout_t layout{666};
    Grid_t dst{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 11};
    Grid_t const src{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 11};
    auto const box = layout.ghostBoxFor(dst);

    {
        ScopeTimer timer{FOR};
        for (std::size_t L = 0; L < FOR; ++L)
        {
            auto d_span       = make_field_box_span(box, dst);
            auto const s_span = make_field_box_span(box, src);

            auto d_slabs    = d_span.begin();
            auto s_slabs    = s_span.begin();
            auto const send = s_span.end();
            for (; s_slabs != send; ++s_slabs, ++d_slabs)
            {
                auto& s_slab = *s_slabs;
                auto& d_slab = *d_slabs;

                auto d_rows      = d_slab.begin();
                auto s_rows      = s_slab.begin();
                auto const srend = s_slab.end();

                for (; s_rows != srend; ++s_rows, ++d_rows)
                {
                    auto& s_row = *s_rows;
                    auto& d_row = *d_rows;

                    for (std::size_t i = 0; i < s_row.size(); ++i)
                        d_row[i] += s_row[i] + L;
                }
            }
        }
    }

    PHARE_DEBUG_DO({
        for (auto const bix : box)
            if (dst(bix) != 11 * 11 + 45)
                return 1;
    })

    return dst.data()[0] != dst.data()[dst.size() - 1];
}

} // namespace PHARE::bench::core

int main()
{
    std::ios::sync_with_stdio(false);
    auto const ret = PHARE::bench::core::test();
    PRINT(ret);
    return ret;
}
