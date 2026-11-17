// Cell deposit (CellDeposit, sorted domain arrays) against the per-particle
// deposit (Interpolator::depositOn): same particles, density and flux equal up
// to the summation order.

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/interpolator/cell_deposit.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/utilities/point/point.hpp"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

namespace PHARE::core
{
// the AMR to local index shift of a patch whose first domain cell is local
// index `ghosts`
template<std::size_t dim>
struct ShiftLayout
{
    static constexpr int ghosts = 5;
    Point<int, dim> AMRToLocal(Point<int, dim> amr) const
    {
        for (std::size_t d = 0; d < dim; ++d)
            amr[d] += ghosts;
        return amr;
    }
};

template<std::size_t dim>
struct Moments
{
    std::array<std::size_t, dim> shape;
    std::vector<double> n, x, y, z;

    explicit Moments(std::size_t const nodes)
    {
        shape.fill(nodes);
        std::size_t size = 1;
        for (auto s : shape)
            size *= s;
        n.assign(size, 0.);
        x.assign(size, 0.);
        y.assign(size, 0.);
        z.assign(size, 0.);
    }

    struct Targets
    {
        double *density, *xFlux, *yFlux, *zFlux;
        std::array<std::size_t, dim> shape;
    };
    Targets targets() { return {n.data(), x.data(), y.data(), z.data(), shape}; }
};

// max |a - b| / max |b| over one moment
inline double relDiff(std::vector<double> const& a, std::vector<double> const& b)
{
    double num = 0, den = 0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        num = std::max(num, std::abs(a[i] - b[i]));
        den = std::max(den, std::abs(b[i]));
    }
    return den > 0 ? num / den : num;
}

enum class ArrayLayout { unsorted, rows, bags };

template<std::size_t dim, std::size_t interp>
void cellDepositMatchesPerParticle(ArrayLayout const layout)
{
#if PHARE_HYB_XSIMD
    constexpr int cells = 6;
    constexpr int ppc   = 21; // not a multiple of the lane count: partial blocks in every row
    std::mt19937_64 gen(1234 + 10 * dim + interp);
    std::uniform_real_distribution<double> unit(0., 1.);
    std::uniform_int_distribution<int> cell(0, cells - 1);

    Box<int, dim> box;
    for (std::size_t d = 0; d < dim; ++d)
    {
        box.lower[d] = 0;
        box.upper[d] = cells - 1;
    }
    ParticleArray<dim> particles(box);
    for (int i = 0; i < ppc * static_cast<int>(box.size()); ++i)
    {
        Particle<dim> p;
        p.weight = .5 + unit(gen);
        for (std::size_t d = 0; d < dim; ++d)
        {
            p.iCell[d] = cell(gen);
            p.delta[d] = unit(gen);
        }
        for (auto& v : p.v)
            v = 2. * unit(gen) - 1.;
        particles.push_back(p);
    }

    typename ParticleArray<dim>::SortScratch scratch;
    auto const keyOf = [&](auto const& p) { return particles.key(p.iCell); };
    if (layout == ArrayLayout::rows)
        particles.layoutRows(keyOf, scratch);
    else if (layout == ArrayLayout::bags)
    {
        std::vector<std::pair<std::size_t, int>> noCrossers;
        particles.sortBags(keyOf, 0, noCrossers, scratch);
        ASSERT_TRUE(particles.bagged());
    }

    std::size_t const nodes = cells + 2 * ShiftLayout<dim>::ghosts + 1;
    Moments<dim> reference(nodes), cellDeposited(nodes);

    ShiftLayout<dim> const shift;
    Interpolator<dim, interp> interpolator;
    auto const refTargets = reference.targets();
    particles.forEachLiveRun([&](std::size_t const first, std::size_t const last) {
        for (auto i = first; i < last; ++i)
            interpolator.depositOn(particles[i], refTargets, shift);
    });

    std::array<int, dim> a2l;
    a2l.fill(ShiftLayout<dim>::ghosts);
    CellDeposit<dim, interp> cellDeposit;
    cellDeposit.depositArray(particles, cellDeposited.targets(), a2l);

    // the b-spline weights of a particle sum to 1: the density sums to the total weight
    double totalWeight = 0, totalDensity = 0;
    particles.forEachLiveRun([&](std::size_t const first, std::size_t const last) {
        for (auto i = first; i < last; ++i)
            totalWeight += particles[i].weight;
    });
    for (auto const n : reference.n)
        totalDensity += n;
    EXPECT_NEAR(totalDensity, totalWeight, 1e-12 * totalWeight);

    double const tol = 1e-14;
    EXPECT_LE(relDiff(cellDeposited.n, reference.n), tol);
    EXPECT_LE(relDiff(cellDeposited.x, reference.x), tol);
    EXPECT_LE(relDiff(cellDeposited.y, reference.y), tol);
    EXPECT_LE(relDiff(cellDeposited.z, reference.z), tol);
#else
    (void)layout;
    GTEST_SKIP() << "cell kernels need the xsimd backend";
#endif
}

template<std::size_t dim, std::size_t interp>
void allLayouts()
{
    for (auto layout : {ArrayLayout::unsorted, ArrayLayout::rows, ArrayLayout::bags})
    {
        SCOPED_TRACE("dim " + std::to_string(dim) + " interp " + std::to_string(interp) + " layout "
                     + std::to_string(static_cast<int>(layout)));
        cellDepositMatchesPerParticle<dim, interp>(layout);
    }
}

TEST(CellDeposit, matchesPerParticleDeposit2D)
{
    allLayouts<2, 1>();
    allLayouts<2, 2>();
    allLayouts<2, 3>();
}

TEST(CellDeposit, matchesPerParticleDeposit3D)
{
    allLayouts<3, 1>();
    allLayouts<3, 2>();
    allLayouts<3, 3>();
}

} // namespace PHARE::core

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
