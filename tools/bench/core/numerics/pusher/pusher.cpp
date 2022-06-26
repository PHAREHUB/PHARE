
#include "benchmark/benchmark.hpp"

#include "phare_core.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"

template<std::size_t dim>
using Field = PHARE::core::Field<PHARE::core::NdArrayVector<dim>,
                                 typename PHARE::core::HybridQuantity::Scalar>;
template<std::size_t dim>
using VecField
    = PHARE::core::VecField<PHARE::core::NdArrayVector<dim>, typename PHARE::core::HybridQuantity>;

template<std::size_t dim>
PHARE::core::Particle<dim> particle()
{
    return {//
            /*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ PHARE::core::ConstArray<int, dim>(35),
            /*.delta  = */ PHARE::core::ConstArray<double, dim>(.01),
            /*.v      = */ {{0, 10., 0}}};
}

template<typename GridLayout, typename Quantity, std::size_t dim = GridLayout::dimension>
Field<dim> field(std::string key, Quantity type, GridLayout const& layout)
{
    Field<dim> feeld{key, type, layout.allocSize(type)};
    std::fill(feeld.begin(), feeld.end(), 1);
    return feeld;
}

template<std::size_t dim, std::size_t interp>
void push(benchmark::State& state)
{
    constexpr std::uint32_t cells = 65;
    constexpr std::uint32_t parts = 1e7;

    using PHARE_Types       = PHARE::core::PHARE_Types<dim, interp>;
    using Interpolator      = PHARE::core::Interpolator<dim, interp>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dim, interp>;
    using Ions_t            = typename PHARE_Types::Ions_t;
    using Electromag_t      = typename PHARE_Types::Electromag_t;
    using GridLayout_t      = typename PHARE_Types::GridLayout_t;
    using ParticleArray     = typename Ions_t::particle_array_type;
    using PartIterator      = typename ParticleArray::iterator;

    using BorisPusher_t = PHARE::core::BorisPusher<dim, PartIterator, Electromag_t, Interpolator,
                                                   BoundaryCondition, GridLayout_t>;

    Interpolator interpolator;
    ParticleArray domainParticles{parts, particle<dim>()};
    ParticleArray tmpDomain{domainParticles.size(), particle<dim>()};

    auto rangeIn  = PHARE::core::makeRange(domainParticles);
    auto rangeOut = PHARE::core::makeRange(tmpDomain);

    auto meshSize = PHARE::core::ConstArray<double, dim>(1.0 / cells);
    auto nCells   = PHARE::core::ConstArray<std::uint32_t, dim>(cells);
    auto origin   = PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)};
    GridLayout_t layout{meshSize, nCells, origin};

    Field<dim> bx = field("Bx", PHARE::core::HybridQuantity::Scalar::Bx, layout);
    Field<dim> by = field("By", PHARE::core::HybridQuantity::Scalar::By, layout);
    Field<dim> bz = field("Bz", PHARE::core::HybridQuantity::Scalar::Bz, layout);

    Field<dim> ex = field("Ex", PHARE::core::HybridQuantity::Scalar::Ex, layout);
    Field<dim> ey = field("Ey", PHARE::core::HybridQuantity::Scalar::Ey, layout);
    Field<dim> ez = field("Ez", PHARE::core::HybridQuantity::Scalar::Ez, layout);

    PHARE::core::Electromag<VecField<dim>> emFields{std::string{"EM"}};
    emFields.B.setBuffer("EM_B_x", &bx);
    emFields.B.setBuffer("EM_B_y", &by);
    emFields.B.setBuffer("EM_B_z", &bz);
    emFields.E.setBuffer("EM_E_x", &ex);
    emFields.E.setBuffer("EM_E_y", &ey);
    emFields.E.setBuffer("EM_E_z", &ez);

    BorisPusher_t pusher;
    pusher.setMeshAndTimeStep(layout.meshSize(), .001);

    while (state.KeepRunning())
    {
        pusher.move(
            /*ParticleRange const&*/ rangeIn, /*ParticleRange&*/ rangeOut,
            /*Electromag const&*/ emFields, /*double mass*/ 1, /*Interpolator&*/ interpolator,
            /*ParticleSelector const&*/ [](auto const& /*part*/) { return true; },
            /*GridLayout const&*/ layout);
    }
}
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/1, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/2, /*interp=*/3)->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/1)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/2)->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(push, /*dim=*/3, /*interp=*/3)->Unit(benchmark::kMicrosecond);

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
