#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_slabs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "benchmark/benchmark.h"

#include <cstddef>
#include <cstdint>


namespace PHARE::core::bench
{
std::uint32_t constexpr static cells = 64;

// for_each_point with its innermost loop kept scalar, as a baseline
template<typename T, std::size_t dim>
void for_each_point_scalar(Box<T, dim> const& box, auto&& fn)
{
    for (auto const& slab : make_box_slabs(box))
        for (auto const& span : slab)
#if defined(__clang__)
#pragma clang loop vectorize(disable) interleave(disable)
#elif defined(__GNUC__)
#pragma GCC novector
#endif
            for (auto const& point : span)
                fn(point);
}

// as GridLayout::evalOnAnyBox, the main caller
template<std::size_t dim>
auto unpacked(auto& fn)
{
    return [&](auto const& point) {
        if constexpr (dim == 1)
            fn(point[0]);
        else if constexpr (dim == 2)
            fn(point[0], point[1]);
        else
            fn(point[0], point[1], point[2]);
    };
}

template<std::size_t dim>
auto make_box()
{
    return Box<std::uint32_t, dim>{ConstArray<std::uint32_t, dim>(1),
                                   ConstArray<std::uint32_t, dim>(cells - 2)};
}

template<std::size_t dim>
auto make_array(double const v)
{
    NdArrayVector<dim> arr{ConstArray<std::uint32_t, dim>(cells)};
    std::fill(arr.data(), arr.data() + arr.size(), v);
    return arr;
}


template<std::size_t dim, bool vectorized>
void for_each_point_axpy(benchmark::State& state)
{
    auto box = make_box<dim>();
    auto dst = make_array<dim>(0);
    auto a   = make_array<dim>(1);
    auto b   = make_array<dim>(2);

    // hide the sizes from the optimizer: real boxes and shapes are only known at runtime
    benchmark::DoNotOptimize(box);
    benchmark::DoNotOptimize(dst);
    benchmark::DoNotOptimize(a);
    benchmark::DoNotOptimize(b);

    auto const& ca = a;
    auto const& cb = b;
    auto kernel    = [&](auto const... ijk) { dst(ijk...) = ca(ijk...) + 2 * cb(ijk...); };
    auto fn        = unpacked<dim>(kernel);

    for (auto _ : state)
    {
        if constexpr (vectorized)
            for_each_point(box, fn);
        else
            for_each_point_scalar(box, fn);
        benchmark::DoNotOptimize(dst.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * box.size());
}

BENCHMARK_TEMPLATE(for_each_point_axpy, 1, /*vectorized*/ false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(for_each_point_axpy, 1, /*vectorized*/ true)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(for_each_point_axpy, 2, /*vectorized*/ false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(for_each_point_axpy, 2, /*vectorized*/ true)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(for_each_point_axpy, 3, /*vectorized*/ false)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(for_each_point_axpy, 3, /*vectorized*/ true)->Unit(benchmark::kMicrosecond);

} // namespace PHARE::core::bench


int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
