// collisions between the children of split particles: siblings copy their parent's velocity and
// share a weight, so their pairs would have nothing to exchange. PHARE's patterns put the children
// of a parent in different fine cells, so siblings hardly meet. On fine cells made of real splits,
// with parents of different weights, the momentum and energy fix must converge.

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/numerics/collisions/collisions.hpp"
#include "amr/data/particles/refine/split.hpp"

#include "tests/core/numerics/collisions/test_collisions_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <random>

using namespace PHARE::core;

namespace
{
template<std::size_t dim, std::size_t interp, std::size_t nbrRefinedParts>
using Splitter = PHARE::amr::Splitter<DimConst<dim>, InterpConst<interp>,
                                      RefinedParticlesConst<nbrRefinedParts>>;

template<typename Splitter_t, std::size_t dim_, std::size_t coarseCells_, std::size_t ppc_>
struct Case
{
    using Splitter                           = Splitter_t;
    static constexpr std::size_t dim         = dim_;
    static constexpr std::size_t coarseCells = coarseCells_;
    static constexpr std::size_t ppc         = ppc_;
};

template<typename Case>
struct SplitCollisionsTest : public ::testing::Test
{
};

using Cases
    = ::testing::Types<Case<Splitter<1, 1, 2>, 1, 32, 64>, Case<Splitter<1, 1, 3>, 1, 32, 64>,
                       Case<Splitter<1, 3, 5>, 1, 32, 64>, Case<Splitter<2, 1, 9>, 2, 8, 32>,
                       Case<Splitter<2, 2, 4>, 2, 8, 32>, Case<Splitter<2, 2, 16>, 2, 8, 16>,
                       Case<Splitter<2, 3, 25>, 2, 8, 16>, Case<Splitter<3, 1, 6>, 3, 4, 16>,
                       Case<Splitter<3, 1, 27>, 3, 4, 16>, Case<Splitter<3, 2, 12>, 3, 4, 16>,
                       Case<Splitter<3, 3, 27>, 3, 4, 16>>;
TYPED_TEST_SUITE(SplitCollisionsTest, Cases);


TYPED_TEST(SplitCollisionsTest, fixConvergesOnSplitParticles)
{
    auto constexpr dim         = TypeParam::dim;
    auto constexpr coarseCells = static_cast<int>(TypeParam::coarseCells);
    typename TypeParam::Splitter const split{};

    // parents with a density varying from cell to cell, so neighbouring families differ in weight
    auto const fineBox = [&]() {
        Box<int, dim> box;
        for (std::size_t d = 0; d < dim; ++d)
            box.lower[d] = 0, box.upper[d] = 2 * coarseCells - 1;
        return box;
    }();
    TestIonsND<dim> ions;
    ions.populations.emplace_back(1., 1., fineBox);
    auto& fine = ions.populations[0].domainParticles();

    std::mt19937_64 rng{7};
    std::uniform_real_distribution<double> uniform{0., 1.};
    std::normal_distribution<double> normal{0., 0.2};
    std::vector<Particle<dim>> children(split.nbRefinedParts);

    Box<int, dim> coarseBox;
    for (std::size_t d = 0; d < dim; ++d)
        coarseBox.lower[d] = 0, coarseBox.upper[d] = coarseCells - 1;
    for (auto const& cell : coarseBox)
    {
        auto const density = 1. + 0.5 * std::sin(2. * M_PI * cell[0] / coarseCells);
        for (std::size_t p = 0; p < TypeParam::ppc; ++p)
        {
            // the parent's position on the fine grid, as particles_data_split.hpp computes it
            // before split()
            Particle<dim> parent;
            parent.weight = density / TypeParam::ppc;
            parent.charge = 1.;
            for (std::size_t d = 0; d < dim; ++d)
            {
                auto const x    = 2. * (cell[d] + uniform(rng));
                parent.iCell[d] = static_cast<int>(std::floor(x));
                parent.delta[d] = x - parent.iCell[d];
            }
            parent.v = {normal(rng), normal(rng), normal(rng)};

            split(parent, children);
            for (auto const& child : children)
                if (isIn(Point{child.iCell}, fineBox))
                    fine.push_back(child);
        }
    }

    auto const totals = [&]() {
        std::array<double, 4> sum{};
        for (auto const& p : fine)
            for (std::size_t i = 0; i < 3; ++i)
            {
                sum[i] += p.weight * p.v[i];
                sum[3] += 0.5 * p.weight * p.v[i] * p.v[i];
            }
        return sum;
    };

    // s at the thermal speed ~0.05, the upper end of the accurate range, where the weighted
    // acceptance leaves the largest energy errors for the fix to absorb
    Collisions<TestIonsND<dim>> collide{{0.03, 1, ScatteringLaw::nanbu, true, 11}};
    auto const before    = totals();
    auto constexpr steps = 20;
    for (std::uint64_t step = 0; step < steps; ++step)
        collide(ions, fineBox, collisions::combine(3, step), 0.01);
    auto const after = totals();

    auto const& counters = collide.counters();
    auto const blocks    = static_cast<double>(fineBox.size() * steps);
    std::cout << dim << "D, " << split.nbRefinedParts << " children, "
              << fine.size() / fineBox.size() << " per fine cell: skipped pairs "
              << static_cast<double>(counters.skipped) / counters.pairs << ", unconverged blocks "
              << counters.unconverged / blocks << "\n";

    // siblings rarely share a fine cell (some 3D patterns, < 0.05% of the pairs)
    EXPECT_LT(static_cast<double>(counters.skipped) / counters.pairs, 1e-3);
    EXPECT_LT(counters.unconverged / blocks, 0.01);
    for (std::size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(after[i], before[i], 1e-12 * before[3]);
    EXPECT_NEAR(after[3] - before[3], counters.residualEnergy, 1e-12 * before[3]);
}

} // namespace


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
