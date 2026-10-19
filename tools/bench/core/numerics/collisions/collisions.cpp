// cost of one collision step per particle, to compare with the pusher benchmark, for the
// particle counts per cell met in PHARE runs, one and two species groups, without and with the
// momentum and energy fix

#include "core/numerics/collisions/collisions.hpp"

#include "tests/core/numerics/collisions/test_collisions_fixtures.hpp"

#include "benchmark/benchmark.h"

#include <random>

using namespace PHARE::core;

void collide(benchmark::State& state)
{
    auto const ppc       = static_cast<std::size_t>(state.range(0));
    auto const nbrGroups = state.range(1);
    auto const exact     = state.range(2) == 1;
    auto constexpr cells = 64;

    std::mt19937_64 rng{1};
    TestIons ions;
    ions.populations.emplace_back(1., 1., cells);
    ions.populations[0].load(ppc, 1., {0., 0., 0.}, {0.2, 0.3, 0.3}, rng, true);
    if (nbrGroups == 2)
    {
        ions.populations.emplace_back(4., 2., cells);
        ions.populations[1].load(ppc, 0.1, {0.1, 0., 0.}, {0.1, 0.1, 0.1}, rng, true);
    }

    Collisions<TestIons> collisions{{1e-3, 1, ScatteringLaw::nanbu, exact, 1}};
    auto const box = Box<int, 1>{{0}, {cells - 1}};

    std::uint64_t step = 0;
    for (auto _ : state)
        collisions(ions, box, collisions::combine(1, step++), 0.01);

    state.SetItemsProcessed(state.iterations() * cells * ppc * nbrGroups);
}
BENCHMARK(collide)->ArgsProduct({{20, 100, 500}, {1, 2}, {0, 1}})->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
