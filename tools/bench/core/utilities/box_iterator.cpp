#include "bench_box.hpp"
#include "core/def.hpp"



namespace PHARE::bench::core
{

int test()
{
    std::size_t constexpr static FOR = 10;

    GridLayout_t layout{666};
    Grid_t dst{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 11};
    Grid_t src{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 11};
    auto const box = layout.ghostBoxFor(dst);

    {
        ScopeTimer timer{FOR};
        for (std::size_t L = 0; L < FOR; ++L)
        {
            auto src_it = box.begin();
            auto dst_it = box.begin();
            for (; dst_it != box.end(); ++src_it, ++dst_it)
                dst(*dst_it) += src(*src_it) + L;
        }
    }

    PHARE_DEBUG_DO({
        for (auto const bix : layout.ghostBoxFor(dst))
            if (dst(bix) != 11 * 11 + 45)
                return 1;
    })

    return dst.data()[0] != dst.data()[dst.size() - 1];
}

} // namespace PHARE::bench::core

int main()
{
    auto const ret = PHARE::bench::core::test();
    PRINT(ret);
    return ret;
}
