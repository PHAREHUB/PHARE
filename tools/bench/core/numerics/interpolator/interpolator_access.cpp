
#include "simulator/simulator_def.hpp"

#include "tools/bench/core/bench.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tools/bench/core/numerics/interpolator/bench_intepolator_access.hpp"

#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <cstddef>
#include <cstdint>


// MeshToParticle and ParticleToMesh timed separately, per version:
//  0 = previous::Interpolator (frozen copy of the per point field access, per field deposit)
//  1 = production Interpolator (stride based field access, fused __restrict deposit)
//  2D/3D only: in 1D the stencil is already a single contiguous line
//  before timing, production is checked to produce the same values as the previous version
//  sorted = particles ordered by iCell (lexicographic, so following the C ordered field
//   memory), otherwise left in their random dispersed order


namespace PHARE::core::bench
{
std::uint32_t constexpr static cells   = 30;
std::uint32_t constexpr static n_parts = 1e6;

template<std::size_t dim, std::size_t interp>
struct AccessSetup
{
    auto static constexpr opts = SimOpts{dim, interp};

    using PHARE_Types   = core::PHARE_Types<opts>;
    using GridLayout_t  = TestGridLayout<typename PHARE_Types::Hybrid::GridLayout_t>;
    using ParticleArray = PHARE_Types::Hybrid::ParticleArray_t;
    using Grid_t        = PHARE_Types::Hybrid::Grid_t;

    AccessSetup(bool const sorted)
    {
        particles.vector().resize(n_parts, particle<dim>());
        disperse(particles, 0, cells - 1, /*seed*/ 1337);

        // non zero weight and random deltas, so the check below compares something
        std::mt19937_64 gen{1337};
        std::uniform_real_distribution<double> distrib{0, 1};
        for (auto& p : particles)
        {
            p.weight = 1;
            for (auto& d : p.delta)
                d = distrib(gen);
        }

        if (sorted)
            std::sort(particles.vector().begin(), particles.vector().end(),
                      [](auto const& a, auto const& b) { return a.iCell < b.iCell; });

        for (std::size_t i = 0; i < 3; ++i)
            for (auto* field : {&(*em).E[i], &(*em).B[i]})
                for (std::size_t j = 0; j < field->size(); ++j)
                    field->data()[j] = std::sin(j + i);
    }

    auto make_density() const
    {
        return Grid_t{"particleDensity", HybridQuantity::Scalar::rho,
                      layout.allocSize(HybridQuantity::Scalar::rho)};
    }

    auto make_flux() const { return UsableVecField<dim>{"F", layout, HybridQuantity::Vector::V}; }

    GridLayout_t layout{cells};
    ParticleArray particles{layout.AMRBox()};
    UsableElectromag<dim> em{layout};
};


bool near(double const a, double const b)
{
    return std::abs(a - b) <= 1e-12 * std::max({1., std::abs(a), std::abs(b)});
}

// out of line, as are the timed passes below: inlined, the check (which runs both versions'
//  ParticleToMesh) changes the register allocation of the timed loop around it
template<std::size_t dim, std::size_t interp>
[[gnu::noinline]] std::string check(AccessSetup<dim, interp>& setup)
{
    Interpolator_t<dim, interp, /*previous*/ 0> previous;
    Interpolator_t<dim, interp, /*production*/ 1> production;

    for (auto& particle : setup.particles)
        if (previous(particle, setup.em, setup.layout)
            != production(particle, setup.em, setup.layout))
        {
            auto const [pE, pB] = previous(particle, setup.em, setup.layout);
            auto const [qE, qB] = production(particle, setup.em, setup.layout);
            for (std::size_t i = 0; i < 3; ++i)
                if (!near(pE[i], qE[i]) or !near(pB[i], qB[i]))
                    return "MeshToParticle mismatch";
        }

    auto p_density = setup.make_density();
    auto q_density = setup.make_density();
    auto p_flux    = setup.make_flux();
    auto q_flux    = setup.make_flux();

    previous(setup.particles, p_density, p_flux, setup.layout);
    production(setup.particles, q_density, q_flux, setup.layout);

    auto const same = [](auto const& a, auto const& b) {
        for (std::size_t i = 0; i < a.size(); ++i)
            if (!near(a.data()[i], b.data()[i]))
                return false;
        return true;
    };

    if (!same(p_density, q_density))
        return "ParticleToMesh density mismatch";
    for (std::size_t i = 0; i < 3; ++i)
        if (!same(p_flux[i], q_flux[i]))
            return "ParticleToMesh flux mismatch";

    return "";
}

template<typename Interpolator_, typename Setup>
[[gnu::noinline]] void mesh_to_particle_pass(Interpolator_& interpolator, Setup& setup)
{
    for (auto& particle : setup.particles)
    {
        auto EB = interpolator(particle, setup.em, setup.layout);
        benchmark::DoNotOptimize(EB);
    }
}

template<typename Interpolator_, typename Setup, typename Density, typename Flux>
[[gnu::noinline]] void particle_to_mesh_pass(Interpolator_& interpolator, Setup& setup,
                                             Density& density, Flux& flux)
{
    interpolator(setup.particles, density, flux, setup.layout);
}


template<std::size_t dim, std::size_t interp, int version, bool sorted>
void mesh_to_particle(benchmark::State& state)
{
    AccessSetup<dim, interp> setup{sorted};
    if (auto const err = check(setup); !err.empty())
        return state.SkipWithError(err.c_str());

    Interpolator_t<dim, interp, version> interpolator;

    for (auto _ : state)
        mesh_to_particle_pass(interpolator, setup);

    state.SetItemsProcessed(state.iterations() * setup.particles.size());
}

template<std::size_t dim, std::size_t interp, int version, bool sorted>
void particle_to_mesh(benchmark::State& state)
{
    AccessSetup<dim, interp> setup{sorted};
    if (auto const err = check(setup); !err.empty())
        return state.SkipWithError(err.c_str());

    Interpolator_t<dim, interp, version> interpolator;
    auto density = setup.make_density();
    auto flux    = setup.make_flux();

    for (auto _ : state)
    {
        particle_to_mesh_pass(interpolator, setup, density, flux);
        benchmark::DoNotOptimize(density.data());
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * setup.particles.size());
}


#define PHARE_ACCESS_BENCH_(fn, dim, interp, sorted)                                               \
    BENCHMARK_TEMPLATE(fn, dim, interp, /*previous*/ 0, sorted)->Unit(benchmark::kMicrosecond);    \
    BENCHMARK_TEMPLATE(fn, dim, interp, /*production*/ 1, sorted)->Unit(benchmark::kMicrosecond);

#define PHARE_ACCESS_BENCH(fn, dim, interp)                                                        \
    PHARE_ACCESS_BENCH_(fn, dim, interp, /*sorted*/ false)                                         \
    PHARE_ACCESS_BENCH_(fn, dim, interp, /*sorted*/ true)

PHARE_ACCESS_BENCH(mesh_to_particle, 2, 1)
PHARE_ACCESS_BENCH(mesh_to_particle, 2, 2)
PHARE_ACCESS_BENCH(mesh_to_particle, 2, 3)
PHARE_ACCESS_BENCH(mesh_to_particle, 3, 1)
PHARE_ACCESS_BENCH(mesh_to_particle, 3, 2)
PHARE_ACCESS_BENCH(mesh_to_particle, 3, 3)

PHARE_ACCESS_BENCH(particle_to_mesh, 2, 1)
PHARE_ACCESS_BENCH(particle_to_mesh, 2, 2)
PHARE_ACCESS_BENCH(particle_to_mesh, 2, 3)
PHARE_ACCESS_BENCH(particle_to_mesh, 3, 1)
PHARE_ACCESS_BENCH(particle_to_mesh, 3, 2)
PHARE_ACCESS_BENCH(particle_to_mesh, 3, 3)

#undef PHARE_ACCESS_BENCH
#undef PHARE_ACCESS_BENCH_

} // namespace PHARE::core::bench


int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
