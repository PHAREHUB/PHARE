#include "bench_box.hpp"
#include "core/def.hpp"



namespace PHARE::bench::core
{

struct Opper
{
    std::size_t L, off;
    double* dst;
    double* src;

    void operator()(auto i)
    {
        auto const idx = off + i;
        dst[idx] += src[idx] + L;
    }
};

int test()
{
    std::size_t constexpr static FOR = 10;
    std::size_t constexpr static N   = 8;

    GridLayout_t layout{666}; //
    Grid_t dst{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 11};
    Grid_t src{"rho", layout, PHARE::core::HybridQuantity::Scalar::rho, 11};

    auto const box   = layout.ghostBoxFor(dst);
    auto const shape = box.shape();
    auto const size  = shape[2];
    auto const mod   = shape[2] & 100;
    auto const times = shape[2] / 100;

    {
        ScopeTimer timer{FOR};
        for (std::size_t L = 0; L < FOR; ++L)
        {
            auto tslabs  = make_field_box_span(box, dst);
            auto rslabs  = make_field_box_span(box, src);
            auto tslabit = tslabs.begin();

            for (auto rslabit = rslabs.begin(); rslabit != rslabs.end(); ++rslabit, ++tslabit)
            {
                auto rrow   = *rslabit;
                auto trow   = *tslabit;
                auto trowit = trow.begin();

                for (auto rrowit = rrow.begin(); rrowit != rrow.end(); ++rrowit, ++trowit)
                {
                    std::size_t offset = 0;

                    for (std::size_t i = 0; i < times; ++i)
                    {
                        Opper opper{L, offset, &(*trowit)[0], &(*rrowit)[0]};
                        PHARE::core::for_N<N>(opper);
                        offset += N;
                    }

                    for (std::size_t i = times * N; i < size; ++i)
                        (*trowit)[i] += (*rrowit)[i] + L;
                }
            }
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
