#include "simulator/simulator_def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field_box.hpp"

#include "tools/bench/core/bench.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

#include <array>
#include <vector>
#include <cstdint>

// FieldBox copies, as in MPI unpack (set_from), pack (append_to) and field to field copies
//  (operate_on_fields), on ghost-slab shaped boxes. A slab thin along z, the fastest varying
//  index, is made of spans only as long as the slab is thick, so per span costs dominate.
// Only the public FieldBox API is used, so this compiles against older trees for comparison.


namespace PHARE::core::bench
{
std::uint32_t constexpr static cells = 40;

// slab thickness and extent, in cells, of the boxes below
std::uint32_t constexpr static thickness = 2;
std::uint32_t constexpr static extent    = 37;

enum class Slab { thin_z, thin_x, full };

template<Slab slab>
auto make_box(std::uint32_t const shift)
{
    auto const lo = 1 + shift;
    auto const hi = lo + extent - 1;
    auto const th = lo + thickness - 1;

    using Array = std::array<std::uint32_t, 3>;
    if constexpr (slab == Slab::thin_z)
        return Box<std::uint32_t, 3>{Array{lo, lo, lo}, Array{hi, hi, th}};
    else if constexpr (slab == Slab::thin_x)
        return Box<std::uint32_t, 3>{Array{lo, lo, lo}, Array{th, hi, hi}};
    else
        return Box<std::uint32_t, 3>{Array{lo, lo, lo}, Array{hi, hi, hi}};
}

struct Fixture
{
    using GridLayout_t = TestGridLayout<typename PHARE_Types<SimOpts{3, 1}>::Hybrid::GridLayout_t>;

    Fixture()
    {
        for (auto* vf : {&src, &dst})
            for (auto& c : *vf)
                std::fill(c.begin(), c.end(), 1.);
    }

    GridLayout_t layout{cells};
    UsableVecField<3> src{"src", layout, HybridQuantity::Vector::E};
    UsableVecField<3> dst{"dst", layout, HybridQuantity::Vector::E};
};


template<Slab slab>
void unpack(benchmark::State& state)
{
    Fixture fx;
    auto const box = make_box<slab>(0);
    FieldBox dst{fx.dst[0], fx.layout, box};
    std::vector<double> const vec(box.size(), 2.);

    for (auto _ : state)
    {
        dst.template set_from<SetEqual<double>>(vec);
        benchmark::DoNotOptimize(fx.dst[0].data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * box.size());
}

template<Slab slab>
void pack(benchmark::State& state)
{
    Fixture fx;
    auto const box = make_box<slab>(0);
    FieldBox src{fx.src[0], fx.layout, box};
    std::vector<double> vec;
    vec.reserve(box.size());

    for (auto _ : state)
    {
        vec.clear();
        src.append_to(vec);
        benchmark::DoNotOptimize(vec.data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * box.size());
}

template<Slab slab>
void copy(benchmark::State& state)
{
    Fixture fx;
    // same shape, shifted by one cell so dst and src cover different regions
    FieldBox dst{fx.dst[0], fx.layout, make_box<slab>(1)};
    FieldBox const src{fx.src[0], fx.layout, make_box<slab>(0)};

    for (auto _ : state)
    {
        operate_on_fields<SetEqual<double>>(dst, src);
        benchmark::DoNotOptimize(fx.dst[0].data());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * src.lcl_box.size());
}

#define PHARE_FIELD_BOX_COPY_CASES(op)                                                             \
    BENCHMARK_TEMPLATE(op, Slab::thin_z)->Name(#op "/thin_z")->Unit(benchmark::kMicrosecond);      \
    BENCHMARK_TEMPLATE(op, Slab::thin_x)->Name(#op "/thin_x")->Unit(benchmark::kMicrosecond);      \
    BENCHMARK_TEMPLATE(op, Slab::full)->Name(#op "/full")->Unit(benchmark::kMicrosecond)

PHARE_FIELD_BOX_COPY_CASES(unpack);
PHARE_FIELD_BOX_COPY_CASES(pack);
PHARE_FIELD_BOX_COPY_CASES(copy);

} // namespace PHARE::core::bench


int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
