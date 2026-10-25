#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/field/field_box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "tools/bench/core/no_vectorize.hpp"

#include "benchmark/benchmark.h"

#include <cstddef>
#include <cstdint>
#include <utility>


// Does operate_on_fields vectorize its per span loop? PlusEquals takes the element-wise
//  path of operate_on_span (SetEqual is a memcpy). vectorized=false is the same loop with
//  vectorization disabled, as a baseline.
// Inspect with llvm-optrec / llvm-mca via mkn, timings alone are only a hint.


namespace PHARE::core::bench
{
std::uint32_t constexpr static cells = 64;

// what operate_on_fields needs from a FieldBox, without a GridLayout
template<typename Array_t, std::size_t dim>
struct ArrayBox
{
    Array_t& field;
    Box<std::uint32_t, dim> lcl_box;
};

// operate_on_fields with operate_on_span's element-wise loop kept scalar
template<typename Operator>
PHARE_BENCH_NO_VECTORIZE void operate_on_fields_scalar(auto& dst, auto const& src)
{
    auto d_span       = make_field_box_slabs(dst.lcl_box, dst.field);
    auto const s_span = make_field_box_slabs(src.lcl_box, src.field);

    auto d_slabs = d_span.begin();
    auto s_slabs = s_span.begin();
    for (; s_slabs != s_span.end(); ++s_slabs, ++d_slabs)
    {
        auto d_spans = d_slabs.begin();
        auto s_spans = s_slabs.begin();
        for (; s_spans != s_slabs.end(); ++s_spans, ++d_spans)
        {
            auto* d_data       = (*d_spans).data();
            auto const* s_data = (*s_spans).data();
            auto const size    = (*d_spans).size();
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#endif
            for (std::size_t i = 0; i < size; ++i)
                Operator{d_data[i]}(s_data[i]);
        }
    }
}

// same shape, shifted by one cell so dst and src cover different regions
template<std::size_t dim>
auto make_box(std::uint32_t const shift)
{
    return Box<std::uint32_t, dim>{ConstArray<std::uint32_t, dim>(1 + shift),
                                   ConstArray<std::uint32_t, dim>(cells - 3 + shift)};
}

template<std::size_t dim>
auto make_array(double const v)
{
    NdArrayVector<dim> arr{ConstArray<std::uint32_t, dim>(cells)};
    std::fill(arr.data(), arr.data() + arr.size(), v);
    return arr;
}


template<std::size_t dim, bool vectorized>
void operate_on_fields_plus_equals(benchmark::State& state)
{
    using Operator = PlusEquals<double>;

    auto dst_arr = make_array<dim>(0);
    auto src_arr = make_array<dim>(1);
    ArrayBox<NdArrayVector<dim>, dim> dst{dst_arr, make_box<dim>(1)};
    ArrayBox<NdArrayVector<dim> const, dim> src{src_arr, make_box<dim>(0)};

    // hide the sizes from the optimizer: real boxes and shapes are only known at runtime
    benchmark::DoNotOptimize(dst_arr);
    benchmark::DoNotOptimize(src_arr);
    benchmark::DoNotOptimize(dst.lcl_box);
    benchmark::DoNotOptimize(src.lcl_box);

    for (auto _ : state)
    {
        if constexpr (vectorized)
            operate_on_fields<Operator>(dst, std::as_const(src));
        else
            operate_on_fields_scalar<Operator>(dst, std::as_const(src));
        benchmark::DoNotOptimize(dst_arr.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * dst.lcl_box.size());
}

BENCHMARK_TEMPLATE(operate_on_fields_plus_equals, 1, /*vectorized*/ false)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(operate_on_fields_plus_equals, 1, /*vectorized*/ true)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(operate_on_fields_plus_equals, 2, /*vectorized*/ false)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(operate_on_fields_plus_equals, 2, /*vectorized*/ true)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(operate_on_fields_plus_equals, 3, /*vectorized*/ false)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(operate_on_fields_plus_equals, 3, /*vectorized*/ true)
    ->Unit(benchmark::kMicrosecond);

} // namespace PHARE::core::bench


int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
