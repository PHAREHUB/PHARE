#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/errors.hpp"
#include "core/data/particles/particle.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <random>
#include <vector>
#include <algorithm>


using namespace PHARE::core;


template<typename InterpOrder>
class Esirkepov1D : public ::testing::Test
{
protected:
    static constexpr std::size_t dim          = 1;
    static constexpr std::size_t interp_order = InterpOrder::value;
    using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;
    using Grid_t       = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using Particle_t   = Particle<dim>;
    using Esirkepov_t  = Esirkepov<dim, interp_order>;

    static constexpr std::size_t nCells = 20;
    static constexpr double dx          = 0.1;
    static constexpr double dt          = 0.01;
    static constexpr double charge      = -1.;
    static constexpr double weight      = 0.01;

    GridLayout_t layout{{{dx}}, {{nCells}}, Point{0.}};
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    Grid_t rho0{"rho0", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho),
                0.};
    Grid_t rho1{"rho1", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho),
                0.};

    std::vector<Particle_t> particles;

    Esirkepov_t esirkepov;

    // random particles anywhere in the domain, |vx| dt < maxCourant dx
    void loadRandom(std::size_t n, double maxCourant)
    {
        std::mt19937_64 gen{1234};
        std::uniform_int_distribution<int> cell{0, static_cast<int>(nCells) - 1};
        std::uniform_real_distribution<double> unit{0., 1.};
        auto const vmax = maxCourant * dx / dt;
        for (std::size_t i = 0; i < n; ++i)
            particles.emplace_back(
                weight, std::array<int, 1>{cell(gen)}, std::array<double, 1>{unit(gen)},
                std::array<double, 3>{vmax * (2 * unit(gen) - 1), vmax * (2 * unit(gen) - 1),
                                      vmax * (2 * unit(gen) - 1)});
    }

    void depositAll()
    {
        esirkepov.depositCharge(particles, rho0, layout, charge, -dt / 2);
        esirkepov.depositCharge(particles, rho1, layout, charge, +dt / 2);
        esirkepov(particles, J.super(), layout, charge, dt);
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
};

using InterpOrders = ::testing::Types<std::integral_constant<std::size_t, 1>,
                                      std::integral_constant<std::size_t, 2>,
                                      std::integral_constant<std::size_t, 3>>;
TYPED_TEST_SUITE(Esirkepov1D, InterpOrders);



// (rho^{n+1/2} - rho^{n-1/2}) / dt + d(Jx)/dx = 0 on every primal node, to round-off.
// dual node i sits between primal nodes i and i+1.
TYPED_TEST(Esirkepov1D, continuityToRoundOff)
{
    this->loadRandom(1000, 0.9);
    this->depositAll();

    auto const& Jx  = this->J(Component::X);
    auto const tol  = 1e-12 * this->maxAbs(this->rho0) / this->dt;
    auto const nRho = this->rho0.size();
    ASSERT_EQ(Jx.size() + 1, nRho);

    for (std::uint32_t i = 1; i + 1 < nRho; ++i)
    {
        auto const drhodt = (this->rho1(i) - this->rho0(i)) / this->dt;
        auto const divJ   = (Jx(i) - Jx(i - 1)) / this->dx;
        EXPECT_NEAR(drhodt + divJ, 0., tol) << " at primal node " << i;
    }
}


TYPED_TEST(Esirkepov1D, restParticleHasNoCurrent)
{
    this->particles.emplace_back(this->weight, std::array<int, 1>{7}, std::array<double, 1>{0.3},
                                 std::array<double, 3>{0., 0., 0.});
    this->depositAll();

    for (auto const& Jc : this->J)
        EXPECT_EQ(this->maxAbs(Jc), 0.);
}


// sum over nodes of Jx is q w vx summed over particles: Esirkepov conserves the first moment
TYPED_TEST(Esirkepov1D, totalCurrentIsParticleMoment)
{
    this->loadRandom(200, 0.9);
    this->depositAll();

    for (std::size_t c = 0; c < 3; ++c)
    {
        double expected = 0;
        for (auto const& p : this->particles)
            expected += this->charge * p.weight * p.v[c];
        EXPECT_NEAR(this->sum(this->J[c]), expected, 1e-12 * std::abs(expected) + 1e-14);
    }
}


TYPED_TEST(Esirkepov1D, fastParticleThrows)
{
    auto const v = 1.01 * this->dx / this->dt;
    this->particles.emplace_back(this->weight, std::array<int, 1>{7}, std::array<double, 1>{0.3},
                                 std::array<double, 3>{-v, 0., 0.});
    EXPECT_THROW(
        this->esirkepov(this->particles, this->J.super(), this->layout, this->charge, this->dt),
        DictionaryException);
}


// particles in the first and last domain cells moving outward at the maximum displacement:
// the whole charge and current must land inside the field arrays
TYPED_TEST(Esirkepov1D, borderParticlesDepositInsideGhostBox)
{
    auto const v    = 0.99 * this->dx / this->dt;
    auto const last = static_cast<int>(this->nCells) - 1;
    this->particles.emplace_back(this->weight, std::array<int, 1>{0}, std::array<double, 1>{0.},
                                 std::array<double, 3>{-v, -v, -v});
    this->particles.emplace_back(this->weight, std::array<int, 1>{last},
                                 std::array<double, 1>{0.999999}, std::array<double, 3>{v, v, v});
    this->depositAll();

    auto const totalCharge = 2 * this->charge * this->weight;
    EXPECT_NEAR(this->sum(this->rho0), totalCharge, 1e-14);
    EXPECT_NEAR(this->sum(this->rho1), totalCharge, 1e-14);
    EXPECT_NEAR(this->sum(this->J(Component::X)), 0., 1e-12);
}


// the charge deposit at a shifted position uses the shapes of the moments: same as the
// Interpolator density deposit of particles moved there, including across cell borders
TYPED_TEST(Esirkepov1D, chargeShapesAreTheMomentShapes)
{
    using Interpolator_t = Interpolator<TestFixture::dim, TestFixture::interp_order>;
    this->loadRandom(300, 0.9);

    auto const shift = 0.5 * this->dt;
    this->esirkepov.depositCharge(this->particles, this->rho0, this->layout, 1., shift);

    auto moved = this->particles;
    for (auto& p : moved)
    {
        auto const pos = p.iCell[0] + p.delta[0] + shift * p.v[0] / this->dx;
        p.iCell[0]     = static_cast<int>(std::floor(pos));
        p.delta[0]     = pos - p.iCell[0];
    }
    UsableVecField<TestFixture::dim> flux{"F", this->layout, HybridQuantity::Vector::V};
    Interpolator_t{}(moved, this->rho1, flux.super(), this->layout);

    for (std::size_t i = 0; i < this->rho0.size(); ++i)
        EXPECT_NEAR(this->rho0(i), this->rho1(i), 1e-14) << " at primal node " << i;
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
