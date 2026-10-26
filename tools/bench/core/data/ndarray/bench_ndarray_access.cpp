// shifted access ( see core/numerics/differential_operators.hpp ) written three ways
//  all three are expected to compile to the same code, and therefore time the same

#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/types.hpp"

#include "benchmark/benchmark.h"

#include <cstddef>
#include <cstdint>
#include <utility>

namespace PHARE::core::bench
{

struct IndexSequence
{
    template<std::size_t d, typename Array, std::size_t dim>
    static double shifted(Array const& array, MeshIndex<dim> const& idx, int const shift)
    {
        return [&]<std::size_t... i>(std::index_sequence<i...>) {
            return array((i == d ? idx[i] + shift : idx[i])...);
        }(std::make_index_sequence<dim>{});
    }
};

struct Fold
{
    template<std::size_t d, typename Array, std::size_t dim>
    static double shifted(Array const& array, MeshIndex<dim> const& idx, int const shift)
    {
        return for_N_fold<dim>(
            [&](auto... i) { return array((i == d ? idx[i] + shift : idx[i])...); });
    }
};

struct MakeArray
{
    template<std::size_t d, typename Array, std::size_t dim>
    static double shifted(Array const& array, MeshIndex<dim> const& idx, int const shift)
    {
        return array(
            for_N_make_array<dim>([&](auto i) { return i == d ? idx[i] + shift : idx[i]; }));
    }
};


template<typename Impl>
void shifted_3d(benchmark::State& state)
{
    auto const n = static_cast<std::uint32_t>(state.range(0));
    NdArrayVector<3> const in{{n, n, n}, 1.};
    NdArrayVector<3> out{n, n, n};

    while (state.KeepRunning())
    {
        for (std::size_t i = 1; i < n - 1; ++i)
            for (std::size_t j = 1; j < n - 1; ++j)
                for (std::size_t k = 1; k < n - 1; ++k)
                {
                    MeshIndex<3> const idx{i, j, k};
                    out(i, j, k) = Impl::template shifted<0>(in, idx, 1)
                                   - Impl::template shifted<0>(in, idx, -1)
                                   + Impl::template shifted<1>(in, idx, 1)
                                   - Impl::template shifted<1>(in, idx, -1)
                                   + Impl::template shifted<2>(in, idx, 1)
                                   - Impl::template shifted<2>(in, idx, -1);
                }
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
}

BENCHMARK_TEMPLATE(shifted_3d, IndexSequence)->Unit(benchmark::kMicrosecond)->Arg(32)->Arg(64);
BENCHMARK_TEMPLATE(shifted_3d, Fold)->Unit(benchmark::kMicrosecond)->Arg(32)->Arg(64);
BENCHMARK_TEMPLATE(shifted_3d, MakeArray)->Unit(benchmark::kMicrosecond)->Arg(32)->Arg(64);

} // namespace PHARE::core::bench


int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
