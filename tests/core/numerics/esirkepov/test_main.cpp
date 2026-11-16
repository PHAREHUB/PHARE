#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/errors.hpp"
#include "core/data/particles/particle.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <optional>
#include <random>
#include <sstream>
#include <vector>
#include <algorithm>


using namespace PHARE::core;


// the cells and mesh sizes differ in every direction, so that an exchanged direction shows,
// and the patch does not start at AMR index 0, so that the local index arithmetic is checked
template<std::size_t dim_, std::size_t interp_order_>
struct EsirkepovSetup
{
    static constexpr std::size_t dim          = dim_;
    static constexpr std::size_t interp_order = interp_order_;

    using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;
    using Grid_t       = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using Particle_t   = Particle<dim>;
    using Esirkepov_t  = Esirkepov<dim, interp_order>;
    using Index_t      = std::array<std::uint32_t, dim>;

    static constexpr std::array<double, 3> allMeshSizes{0.1, 0.13, 0.07};
    static constexpr std::array<std::uint32_t, 3> allCells{20, 17, 13};
    static constexpr std::array<int, 3> allLower{37, -12, 5};

    static constexpr double dt     = 0.01;
    static constexpr double charge = -1.;
    static constexpr double weight = 0.01;

    static auto sub(auto const& all)
    {
        std::array<typename std::decay_t<decltype(all)>::value_type, dim> a;
        std::copy_n(all.begin(), dim, a.begin());
        return a;
    }

    static GridLayout_t makeLayout()
    {
        auto const cells = sub(allCells);
        auto const lower = sub(allLower);
        std::array<int, dim> upper;
        for (std::size_t d = 0; d < dim; ++d)
            upper[d] = lower[d] + static_cast<int>(cells[d]) - 1;
        Point<double, dim> origin;
        for (std::size_t d = 0; d < dim; ++d)
            origin[d] = lower[d] * allMeshSizes[d];
        return GridLayout_t{sub(allMeshSizes), cells, origin, Box<int, dim>{lower, upper}};
    }

    std::array<double, dim> const dl = sub(allMeshSizes);
    GridLayout_t layout              = makeLayout();
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    Grid_t rho0{"rho0", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho),
                0.};
    Grid_t rho1{"rho1", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho),
                0.};

    std::vector<Particle_t> particles;

    Esirkepov_t esirkepov;

    // random particles anywhere in the domain, |v_d| dt < maxCourant dl_d in every direction
    void loadRandom(std::size_t n, double maxCourant, std::uint64_t seed = 1234)
    {
        std::mt19937_64 gen{seed};
        std::uniform_real_distribution<double> unit{0., 1.};
        for (std::size_t i = 0; i < n; ++i)
        {
            std::array<int, dim> iCell;
            std::array<double, dim> delta;
            std::array<double, 3> v;
            for (std::size_t d = 0; d < dim; ++d)
            {
                std::uniform_int_distribution<int> cell{allLower[d],
                                                        allLower[d] + int(allCells[d]) - 1};
                iCell[d] = cell(gen);
                delta[d] = unit(gen);
            }
            for (std::size_t c = 0; c < 3; ++c)
                v[c] = maxCourant * allMeshSizes[c] / dt * (2 * unit(gen) - 1);
            particles.emplace_back(weight, iCell, delta, v);
        }
    }

    void reset()
    {
        for (auto& Jc : J)
            Jc.zero();
        rho0.zero();
        rho1.zero();
    }

    // rho at the two ends of the segment, and the current over it
    void depositAll(double const back = 1., double const front = 1.)
    {
        esirkepov.depositCharge(particles, rho0, layout, charge, -back * dt / 2);
        esirkepov.depositCharge(particles, rho1, layout, charge, +front * dt / 2);
        esirkepov(particles, J.super(), layout, charge, dt, back, front);
    }

    static double sum(auto const& field)
    {
        double s = 0;
        for (std::size_t i = 0; i < field.size(); ++i)
            s += field.data()[i];
        return s;
    }

    static double maxAbs(auto const& field)
    {
        double m = 0;
        for (std::size_t i = 0; i < field.size(); ++i)
            m = std::max(m, std::abs(field.data()[i]));
        return m;
    }

    // calls fn(index) on every node of an array of the given shape
    static void forEachNode(Index_t const& shape, auto&& fn)
    {
        Index_t idx{};
        while (true)
        {
            fn(idx);
            std::size_t d = dim;
            while (d-- > 0)
            {
                if (++idx[d] < shape[d])
                    break;
                idx[d] = 0;
            }
            if (d > dim) // wrapped around in the first direction
                return;
        }
    }

    // largest |d rho / dt + div J| on the primal nodes where every dual neighbour exists,
    // with dual node i between primal nodes i and i+1 in its direction,
    // in the box of the given extent starting at lower (default: everywhere)
    double maxContinuityResidual(double const duration, Index_t const lower = {},
                                 std::optional<Index_t> const extent = std::nullopt) const
    {
        double residual = 0;
        forEachNode(extent ? *extent : rho0.shape(), [&](Index_t idx) {
            for (std::size_t d = 0; d < dim; ++d)
                idx[d] += lower[d];
            double divJ = 0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                auto const& Jd = J[d];
                if (idx[d] == 0 or idx[d] >= Jd.shape()[d])
                    return;
                auto left = idx;
                --left[d];
                divJ += (Jd(idx) - Jd(left)) / dl[d];
            }
            residual = std::max(residual, std::abs((rho1(idx) - rho0(idx)) / duration + divJ));
        });
        return residual;
    }
};


template<std::size_t dim_, std::size_t interp_order_>
struct Case
{
    static constexpr std::size_t dim          = dim_;
    static constexpr std::size_t interp_order = interp_order_;
};

template<typename Case_t>
class EsirkepovTest : public ::testing::Test,
                      public EsirkepovSetup<Case_t::dim, Case_t::interp_order>
{
};

using Cases = ::testing::Types<Case<1, 1>, Case<1, 2>, Case<1, 3>, Case<2, 1>, Case<2, 2>,
                               Case<2, 3>, Case<3, 1>, Case<3, 2>, Case<3, 3>>;
TYPED_TEST_SUITE(EsirkepovTest, Cases);



TYPED_TEST(EsirkepovTest, currentIsYeeCentered)
{
    for (std::size_t c = 0; c < 3; ++c)
        for (std::size_t d = 0; d < TestFixture::dim; ++d)
            EXPECT_EQ(this->J[c].shape()[d] + (c == d ? 1 : 0), this->rho0.shape()[d])
                << " component " << c << " direction " << d;
}


// (rho^{n+1/2} - rho^{n-1/2}) / dt + div J = 0 on every primal node, to round-off,
// for the leapfrog segment (back = front = 1) and for the start-up one (back = 1, front = 0)
TYPED_TEST(EsirkepovTest, continuityToRoundOff)
{
    this->loadRandom(1000, 0.9);

    for (auto const [back, front] : {std::pair{1., 1.}, std::pair{1., 0.}})
    {
        this->reset();
        this->depositAll(back, front);

        auto const duration = (back + front) * this->dt / 2;
        auto const tol      = 1e-12 * this->maxAbs(this->rho0) / duration;
        EXPECT_LT(this->maxContinuityResidual(duration), tol)
            << " back " << back << " front " << front;
    }
}


// one particle at a time, starting near the cell borders and the cell middles (where the start
// index changes at order 2), moving across zero, one, two or three of them in every direction
TYPED_TEST(EsirkepovTest, continuityOfEachCrossingParticle)
{
    constexpr auto dim = TestFixture::dim;
    std::array<double, 4> const deltas{0.02, 0.48, 0.52, 0.98};
    std::array<double, 3> const courants{-0.9, 0., 0.9};

    std::array<std::uint32_t, dim> choices;
    choices.fill(deltas.size() * courants.size());

    TestFixture::forEachNode(choices, [&](auto const& choice) {
        std::array<int, dim> iCell;
        std::array<double, dim> delta;
        std::array<double, 3> v{0.3, -0.2, 0.1};
        for (std::size_t d = 0; d < dim; ++d)
        {
            iCell[d] = this->allLower[d] + 5;
            delta[d] = deltas[choice[d] % deltas.size()];
            v[d]     = courants[choice[d] / deltas.size()] * this->dl[d] / this->dt;
        }
        this->particles = {typename TestFixture::Particle_t{this->weight, iCell, delta, v}};
        this->reset();
        this->depositAll();

        // the deposit spans at most 5 nodes around the particle cell
        auto lower = this->layout.AMRToLocal(Point{iCell}).toArray();
        for (auto& l : lower)
            l -= 3;
        typename TestFixture::Index_t extent;
        extent.fill(8);

        auto const tol = 1e-12 * this->maxAbs(this->rho0) / this->dt;
        EXPECT_LT(this->maxContinuityResidual(this->dt, lower, extent), tol)
            << " delta " << delta[0] << " v " << v[0] << " " << v[1] << " " << v[2];
    });
}


TYPED_TEST(EsirkepovTest, restParticleHasNoCurrent)
{
    std::array<int, TestFixture::dim> iCell;
    std::array<double, TestFixture::dim> delta;
    for (std::size_t d = 0; d < TestFixture::dim; ++d)
    {
        iCell[d] = this->allLower[d] + 7;
        delta[d] = 0.3;
    }
    this->particles.emplace_back(this->weight, iCell, delta, std::array<double, 3>{0., 0., 0.});
    this->depositAll();

    for (auto const& Jc : this->J)
        EXPECT_EQ(this->maxAbs(Jc), 0.);
}


// sum over nodes of J is q w v summed over particles: Esirkepov conserves the first moment
TYPED_TEST(EsirkepovTest, totalCurrentIsParticleMoment)
{
    this->loadRandom(200, 0.9);
    this->depositAll();

    for (std::size_t c = 0; c < 3; ++c)
    {
        double expected = 0;
        for (auto const& p : this->particles)
            expected += this->charge * p.weight * p.v[c];
        EXPECT_NEAR(this->sum(this->J[c]), expected, 1e-12 * std::abs(expected) + 1e-14)
            << " component " << c;
    }
}


// a displacement of one cell or more in any direction of the grid throws, but not along an
// invariant direction
TYPED_TEST(EsirkepovTest, fastParticleThrows)
{
    constexpr auto dim = TestFixture::dim;
    std::array<int, dim> iCell;
    std::array<double, dim> delta;
    for (std::size_t d = 0; d < dim; ++d)
    {
        iCell[d] = this->allLower[d] + 7;
        delta[d] = 0.3;
    }

    for (std::size_t c = 0; c < 3; ++c)
    {
        std::array<double, 3> v{0., 0., 0.};
        v[c]            = -1.01 * this->allMeshSizes[c] / this->dt;
        this->particles = {typename TestFixture::Particle_t{this->weight, iCell, delta, v}};
        if (c < dim)
            EXPECT_THROW(this->esirkepov(this->particles, this->J.super(), this->layout,
                                         this->charge, this->dt),
                         DictionaryException)
                << " direction " << c;
        else
            EXPECT_NO_THROW(this->esirkepov(this->particles, this->J.super(), this->layout,
                                            this->charge, this->dt))
                << " direction " << c;
    }
}


// particles in the corner cells of the domain moving outward at the maximum displacement:
// the whole charge and current must land inside the field arrays (checked by the asserts of
// a debug build)
TYPED_TEST(EsirkepovTest, cornerParticlesDepositInsideGhostBox)
{
    constexpr auto dim = TestFixture::dim;
    std::array<std::uint32_t, dim> corners;
    corners.fill(2);

    TestFixture::forEachNode(corners, [&](auto const& corner) {
        std::array<int, dim> iCell;
        std::array<double, dim> delta;
        std::array<double, 3> v{0.5, -0.5, 0.5};
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const upper = corner[d] == 1;
            iCell[d]         = this->allLower[d] + (upper ? int(this->allCells[d]) - 1 : 0);
            delta[d]         = upper ? 0.999999 : 0.;
            v[d]             = (upper ? 0.99 : -0.99) * this->dl[d] / this->dt;
        }
        this->particles.emplace_back(this->weight, iCell, delta, v);
    });
    this->depositAll();

    auto const totalCharge = this->particles.size() * this->charge * this->weight;
    EXPECT_NEAR(this->sum(this->rho0), totalCharge, 1e-14);
    EXPECT_NEAR(this->sum(this->rho1), totalCharge, 1e-14);
    for (std::size_t c = 0; c < 3; ++c)
    {
        double expected = 0;
        for (auto const& p : this->particles)
            expected += this->charge * p.weight * p.v[c];
        EXPECT_NEAR(this->sum(this->J[c]), expected, 1e-12 * std::abs(expected) + 1e-14)
            << " component " << c;
    }
}


// the charge deposit at a shifted position uses the shapes of the moments: same as the
// Interpolator density deposit of particles moved there, including across cell borders
TYPED_TEST(EsirkepovTest, chargeShapesAreTheMomentShapes)
{
    constexpr auto dim   = TestFixture::dim;
    using Interpolator_t = Interpolator<dim, TestFixture::interp_order>;
    this->loadRandom(300, 0.9);

    auto const shift = 0.5 * this->dt;
    this->esirkepov.depositCharge(this->particles, this->rho0, this->layout, 1., shift);

    auto moved = this->particles;
    for (auto& p : moved)
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const pos = p.iCell[d] + p.delta[d] + shift * p.v[d] / this->dl[d];
            p.iCell[d]     = static_cast<int>(std::floor(pos));
            p.delta[d]     = pos - p.iCell[d];
        }
    UsableVecField<dim> flux{"F", this->layout, HybridQuantity::Vector::V};
    Interpolator_t{}(moved, this->rho1, flux.super(), this->layout);

    TestFixture::forEachNode(this->rho0.shape(), [&](auto const& idx) {
        EXPECT_NEAR(this->rho0(idx), this->rho1(idx), 1e-14) << " at primal node " << idx[0];
    });
}



// b-spline of order p at a distance d of the particle, in cells: the shapes of the moments
double bspline(std::size_t const order, double const distance)
{
    auto const d = std::abs(distance);
    if (order == 1)
        return std::max(0., 1. - d);
    if (order == 2)
        return d < 0.5 ? 0.75 - d * d : (d < 1.5 ? 0.5 * (1.5 - d) * (1.5 - d) : 0.);
    return d < 1. ? 2. / 3. - d * d + 0.5 * d * d * d
                  : (d < 2. ? (2. - d) * (2. - d) * (2. - d) / 6. : 0.);
}


// Esirkepov 2001 written as published, independently of Esirkepov<dim, order>: shapes S0 and S1
// at the two ends of the segment on every primal node of the patch, the W decomposition
//   1D  W1 = S1x - S0x
//   2D  W1 = (S1x - S0x)(S0y + S1y) / 2, W2 = (S1y - S0y)(S0x + S1x) / 2
//   3D  W1 = S1x S1y S1z / 3 - S0x S1y S1z / 3 + S1x S0y S1z / 6 - S0x S0y S1z / 6
//          + S1x S1y S0z / 6 - S0x S1y S0z / 6 + S1x S0y S0z / 3 - S0x S0y S0z / 3,
//       W2, W3 by circular permutation,
// J_d(i + 1/2) - J_d(i - 1/2) = -q w dl_d / duration W_d(i), summed from the first node of the
// patch, and the current along an invariant direction
//   1D  J = q w v (S0x + S1x) / 2
//   2D  J = q w v (S1x S1y / 3 + S1x S0y / 6 + S0x S1y / 6 + S0x S0y / 3)
template<typename Setup>
void publishedEsirkepov(Setup const& s, UsableVecField<Setup::dim>& Jref, double const back,
                        double const front)
{
    constexpr auto dim  = Setup::dim;
    using Index_t       = Setup::Index_t;
    using Grid_t        = Setup::Grid_t;
    auto const& layout  = s.layout;
    auto const duration = (back + front) * s.dt / 2;
    auto const primal   = layout.allocSize(HybridQuantity::Scalar::rho);

    // W_d and the invariant direction currents, on the primal nodes
    std::array<Grid_t, 3> W{Grid_t{"W1", HybridQuantity::Scalar::rho, primal, 0.},
                            Grid_t{"W2", HybridQuantity::Scalar::rho, primal, 0.},
                            Grid_t{"W3", HybridQuantity::Scalar::rho, primal, 0.}};

    for (auto const& p : s.particles)
    {
        std::array<double, dim> x0, x1;
        for (std::size_t d = 0; d < dim; ++d)
        {
            auto const local = p.iCell[d] - layout.AMRBox().lower[d]
                               + static_cast<int>(layout.physicalStartIndex(
                                   QtyCentering::primal, static_cast<Direction>(d)));
            x0[d] = local + p.delta[d] - back * s.dt / 2 * p.v[d] / s.dl[d];
            x1[d] = local + p.delta[d] + front * s.dt / 2 * p.v[d] / s.dl[d];
        }
        auto const qw = s.charge * p.weight;

        Setup::forEachNode(primal, [&](Index_t const& idx) {
            std::array<double, 3> S0{1., 1., 1.}, S1{1., 1., 1.};
            for (std::size_t d = 0; d < dim; ++d)
            {
                S0[d] = bspline(Setup::interp_order, idx[d] - x0[d]);
                S1[d] = bspline(Setup::interp_order, idx[d] - x1[d]);
            }
            auto const& [S0x, S0y, S0z] = S0;
            auto const& [S1x, S1y, S1z] = S1;

            std::array<double, 3> w{};
            if constexpr (dim == 1)
            {
                w[0] = S1x - S0x;
                w[1] = p.v[1] * (S0x + S1x) / 2;
                w[2] = p.v[2] * (S0x + S1x) / 2;
            }
            else if constexpr (dim == 2)
            {
                w[0] = (S1x - S0x) * (S0y + S1y) / 2;
                w[1] = (S1y - S0y) * (S0x + S1x) / 2;
                w[2] = p.v[2] * (S1x * S1y / 3 + S1x * S0y / 6 + S0x * S1y / 6 + S0x * S0y / 3);
            }
            else
            {
                w[0] = S1x * S1y * S1z / 3 - S0x * S1y * S1z / 3 + S1x * S0y * S1z / 6
                       - S0x * S0y * S1z / 6 + S1x * S1y * S0z / 6 - S0x * S1y * S0z / 6
                       + S1x * S0y * S0z / 3 - S0x * S0y * S0z / 3;
                w[1] = S1y * S1z * S1x / 3 - S0y * S1z * S1x / 3 + S1y * S0z * S1x / 6
                       - S0y * S0z * S1x / 6 + S1y * S1z * S0x / 6 - S0y * S1z * S0x / 6
                       + S1y * S0z * S0x / 3 - S0y * S0z * S0x / 3;
                w[2] = S1z * S1x * S1y / 3 - S0z * S1x * S1y / 3 + S1z * S0x * S1y / 6
                       - S0z * S0x * S1y / 6 + S1z * S1x * S0y / 6 - S0z * S1x * S0y / 6
                       + S1z * S0x * S0y / 3 - S0z * S0x * S0y / 3;
            }
            for (std::size_t c = 0; c < 3; ++c)
                W[c](idx) += c < dim ? -qw * s.dl[c] / duration * w[c] : qw * w[c];
        });
    }

    for (std::size_t c = 0; c < 3; ++c)
    {
        auto& Jc = Jref[c];
        Setup::forEachNode(Jc.shape(), [&](Index_t const& idx) {
            if (c >= dim)
            {
                Jc(idx) = W[c](idx);
                return;
            }
            double running = 0.;
            auto node      = idx;
            for (node[c] = 0; node[c] <= idx[c]; ++node[c])
                running += W[c](node);
            Jc(idx) = running;
        });
    }
}


// the rotational part of the current is fixed by neither continuity nor the moments: compare
// with the published formulas node by node
TYPED_TEST(EsirkepovTest, sameAsPublishedEsirkepov)
{
    this->loadRandom(TestFixture::dim == 3 ? 50 : 200, 0.9);

    for (auto const [back, front] : {std::pair{1., 1.}, std::pair{1., 0.}})
    {
        this->reset();
        this->esirkepov(this->particles, this->J.super(), this->layout, this->charge, this->dt,
                        back, front);
        UsableVecField<TestFixture::dim> Jref{"Jref", this->layout, HybridQuantity::Vector::J};
        publishedEsirkepov(*this, Jref, back, front);

        for (std::size_t c = 0; c < 3; ++c)
        {
            double maxDiff = 0.;
            TestFixture::forEachNode(this->J[c].shape(), [&](auto const& idx) {
                maxDiff = std::max(maxDiff, std::abs(this->J[c](idx) - Jref[c](idx)));
            });
            auto const relative = maxDiff / this->maxAbs(Jref[c]);
            ::testing::Test::RecordProperty(
                "J" + std::to_string(c) + "_back" + std::to_string(int(back)) + "_front"
                    + std::to_string(int(front)),
                (std::ostringstream{} << std::scientific << relative).str());
            EXPECT_LT(relative, 1e-14)
                << " component " << c << " back " << back << " front " << front;
        }
    }
}



// The current along a direction the grid does not resolve is the current of the next dimension
// summed along that direction: the 1D Jy and Jz are the 2D ones summed over y, the 2D Jz is the
// 3D one summed over z. The current along the resolved directions is fixed by the continuity
// equation, this checks the weights of the others.
template<typename InterpOrder>
class EsirkepovProjection : public ::testing::Test
{
protected:
    static constexpr std::size_t interp_order = InterpOrder::value;

    // the higher dimensional current summed along its last direction, against the lower one
    template<std::size_t dim>
    static void checkProjection(std::size_t const component)
    {
        EsirkepovSetup<dim, interp_order> low;
        EsirkepovSetup<dim + 1, interp_order> high;

        high.loadRandom(200, 0.9);
        for (auto const& p : high.particles)
        {
            std::array<int, dim> iCell;
            std::array<double, dim> delta;
            std::copy_n(p.iCell.begin(), dim, iCell.begin());
            std::copy_n(p.delta.begin(), dim, delta.begin());
            low.particles.emplace_back(p.weight, iCell, delta, p.v);
        }
        low.depositAll();
        high.depositAll();

        auto const& Jlow  = low.J[component];
        auto const& Jhigh = high.J[component];
        auto const tol    = 1e-13 * low.maxAbs(Jlow);

        EsirkepovSetup<dim, interp_order>::forEachNode(Jlow.shape(), [&](auto const& idx) {
            std::array<std::uint32_t, dim + 1> hidx;
            std::copy_n(idx.begin(), dim, hidx.begin());
            double summed = 0;
            for (hidx[dim] = 0; hidx[dim] < Jhigh.shape()[dim]; ++hidx[dim])
                summed += Jhigh(hidx);
            EXPECT_NEAR(Jlow(idx), summed, tol)
                << " component " << component << " at " << idx[0] << " in " << dim << "D";
        });
    }
};

using InterpOrders = ::testing::Types<std::integral_constant<std::size_t, 1>,
                                      std::integral_constant<std::size_t, 2>,
                                      std::integral_constant<std::size_t, 3>>;
TYPED_TEST_SUITE(EsirkepovProjection, InterpOrders);


TYPED_TEST(EsirkepovProjection, invariantDirectionCurrentIsSummedCurrent)
{
    TestFixture::template checkProjection<1>(1); // 1D Jy
    TestFixture::template checkProjection<1>(2); // 1D Jz
    TestFixture::template checkProjection<2>(2); // 2D Jz
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
