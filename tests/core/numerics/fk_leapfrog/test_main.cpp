// The fully kinetic leapfrog on one grid, without communications:
// particles (x, v) and B at integer times, E at half integer times.
//   1. J^n by Esirkepov over x^n -+ (dt/2) v^n
//   2. E^{n+1/2} = E^{n-1/2} + dt c^2 (curl B^n - J^n)
//   3. B^{n+1/2} = B^n - (dt/2) curl E^{n+1/2}
//   4. Boris push with E^{n+1/2}, B^{n+1/2}
//   5. B^{n+1} = B^n - dt curl E^{n+1/2}
// started from E^{-1/2} = E^0 - (dt/2) c^2 (curl B^0 - J_back), J_back over x^0 - (dt/2) v^0 ->
// x^0. Since x^{n+1} - (dt/2) v^{n+1} is the half push position x^{n+1/2}, the Gauss residual div
// E^{n+1/2} - c^2 rho^{n+1/2} must stay at its initial value to round-off.

#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/numerics/ampere/ampere_maxwell.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"
#include "core/numerics/faraday/faraday.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/utilities/range/range.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <numbers>
#include <algorithm>


using namespace PHARE::core;


struct EM
{
    VecField_t<1>& E;
    VecField_t<1>& B;
};


template<typename InterpOrder>
class FKLeapfrog1D : public ::testing::Test
{
protected:
    static constexpr std::size_t dim          = 1;
    static constexpr std::size_t interp_order = InterpOrder::value;
    using GridLayout_t    = PHARE_Types<PHARE::SimOpts{dim, interp_order}>::Hybrid::GridLayout_t;
    using Grid_t          = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;
    using ParticleArray_t = ParticleArray<dim>;
    using Interpolator_t  = Interpolator<dim, interp_order>;
    using Pusher_t        = BorisPusher<dim, IndexRange<ParticleArray_t>, EM, Interpolator_t,
                                        BoundaryCondition<dim, interp_order>, GridLayout_t>;

    static constexpr std::size_t nCells = 64;
    static constexpr double dx          = 0.1;
    static constexpr double c2          = 100.;
    static constexpr double dt          = 0.5 * dx / 10.; // c dt / dx = 0.5
    static constexpr std::size_t nSteps = 50;

    struct Species
    {
        double charge, mass, vmax;
        ParticleArray_t particles;
    };

    GridLayout_t layout{{{dx}}, {{nCells}}, Point{0.}};
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
    UsableVecField<dim> J{"J", layout, HybridQuantity::Vector::J};
    UsableVecField<dim> Bhalf{"Bhalf", layout, HybridQuantity::Vector::B};
    Grid_t rho{"rho", HybridQuantity::Scalar::rho, layout.allocSize(HybridQuantity::Scalar::rho),
               0.};

    std::vector<Species> species;
    Esirkepov<dim, interp_order> esirkepov;
    Interpolator_t interpolator;
    Pusher_t pusher;

    FKLeapfrog1D()
    {
        for (auto& f : E)
            fill_(f, [](double) { return 0.; });
        for (auto& f : B)
            fill_(f, [](double) { return 0.; });
        pusher.setMeshAndTimeStep({dx}, dt);
    }

    // warm ions and electrons in the central cells, with a transverse B
    void loadWarmPlasma()
    {
        auto const box = grow(layout.AMRBox(), 2);
        // |v| dt < 0.3 dx and particles loaded in the central cells: they never get closer
        // than 8 cells to the box edges in nSteps
        species.push_back({+1., 1., 0.5, ParticleArray_t{box}});
        species.push_back({-1., 0.04, 3., ParticleArray_t{box}});

        std::mt19937_64 gen{42};
        std::uniform_int_distribution<int> cell{24, 39};
        std::uniform_real_distribution<double> unit{0., 1.};
        for (auto& s : species)
            for (std::size_t i = 0; i < 300; ++i)
                s.particles.emplace_back(
                    Particle<dim>{0.01,
                                  {cell(gen)},
                                  {unit(gen)},
                                  {s.vmax * (2 * unit(gen) - 1), s.vmax * (2 * unit(gen) - 1),
                                   s.vmax * (2 * unit(gen) - 1)}});

        auto const k = 2 * std::numbers::pi / (nCells * dx);
        fill_(B(Component::X), [](double) { return 1.; });
        fill_(B(Component::Y), [&](double x) { return 0.1 * std::sin(k * x); });
        fill_(B(Component::Z), [&](double x) { return 0.1 * std::cos(k * x); });
    }

    // E^{-1/2} from the backward half segment
    void startUp()
    {
        depositCurrent(1., 0.);
        AmpereMaxwell{layout, c2}(B.super(), J.super(), E.super(), -dt / 2);
    }

    void fill_(auto& f, auto&& fn)
    {
        auto const start  = layout.physicalStartIndex(f, Direction::X);
        auto const isDual = layout.centering(f.physicalQuantity())[0] == QtyCentering::dual;
        layout.evalOnGhostBox(f, [&](auto const& i) {
            f(i) = fn((static_cast<double>(i) - start + (isDual ? 0.5 : 0.)) * dx);
        });
    }

    void zero_(auto& f)
    {
        fill_(f, [](double) { return 0.; });
    }

    void depositCurrent(double back, double front)
    {
        for (auto& f : J)
            zero_(f);
        for (auto const& s : species)
            esirkepov(s.particles, J.super(), layout, s.charge, dt, back, front);
    }

    // div E - c^2 rho, rho deposited at x + timeShift v, on primal nodes whose two neighbor
    // dual nodes are physical
    std::vector<double> gaussResidual(double const timeShift = -dt / 2)
    {
        zero_(rho);
        for (auto const& s : species)
            esirkepov.depositCharge(s.particles, rho, layout, s.charge, timeShift);

        auto const& Ex = E(Component::X);
        std::vector<double> residual;
        auto const first = layout.physicalStartIndex(QtyCentering::primal, Direction::X) + 1;
        auto const last  = layout.physicalEndIndex(QtyCentering::primal, Direction::X) - 1;
        for (auto i = first; i <= last; ++i)
            residual.push_back((Ex(i) - Ex(i - 1)) / dx - c2 * rho(i));
        return residual;
    }

    double maxChargeTerm()
    {
        double m = 0;
        for (std::size_t i = 0; i < rho.size(); ++i)
            m = std::max(m, c2 * std::abs(rho.data()[i]));
        return m;
    }

    void step()
    {
        depositCurrent(1., 1.);
        AmpereMaxwell{layout, c2}(B.super(), J.super(), E.super(), dt);
        Faraday{layout}(B.super(), E.super(), Bhalf.super(), dt / 2);

        for (auto& s : species)
        {
            auto range = makeIndexRange(s.particles);
            auto noop  = [](auto& r) { return r; };
            pusher.move(range, range, EM{E.super(), Bhalf.super()}, s.charge, s.mass, interpolator,
                        layout, noop, noop);
        }

        Faraday{layout}(B.super(), E.super(), B.super(), dt);
    }
};

using InterpOrders = ::testing::Types<std::integral_constant<std::size_t, 1>,
                                      std::integral_constant<std::size_t, 2>,
                                      std::integral_constant<std::size_t, 3>>;
TYPED_TEST_SUITE(FKLeapfrog1D, InterpOrders);



TYPED_TEST(FKLeapfrog1D, gaussResidualStaysAtItsInitialValue)
{
    this->loadWarmPlasma();

    // residual of E^0 with the charge at x^0
    auto const atZero = this->gaussResidual(0.);

    this->startUp();

    auto const initial = this->gaussResidual();
    auto const tol     = 1e-10 * this->maxChargeTerm();
    ASSERT_GT(tol, 0.);

    for (std::size_t i = 0; i < initial.size(); ++i)
        ASSERT_NEAR(initial[i], atZero[i], tol) << "start-up residual differs at node " << i;

    double maxDrift = 0;

    for (std::size_t n = 0; n < this->nSteps; ++n)
    {
        this->step();
        auto const residual = this->gaussResidual();

        double drift = 0;
        for (std::size_t i = 0; i < residual.size(); ++i)
            drift = std::max(drift, std::abs(residual[i] - initial[i]));
        ASSERT_LT(drift, tol) << "at step " << n;
        maxDrift = std::max(maxDrift, drift);
    }
    std::cout << "max Gauss residual drift / max c2 rho = " << maxDrift / (tol / 1e-10) << "\n";
}



// cold electrons on immobile ions, paired loading (no initial charge), with a velocity profile
// vanishing at the slab edges: in 1D every electron sheet oscillates at the plasma frequency,
// here sqrt(c^2 n q^2 / m), at the leapfrog frequency (2/dt) asin(omega_pe dt / 2)
TYPED_TEST(FKLeapfrog1D, coldPlasmaOscillatesAtThePlasmaFrequency)
{
    auto const box   = grow(this->layout.AMRBox(), 2);
    auto const dx    = this->dx;
    auto const dt    = this->dt;
    double const me  = 0.04;
    int const first  = 16;
    int const last   = 47;
    int const ppc    = 20;
    double const v0  = 0.1;
    auto const slab  = (last + 1 - first) * dx;
    auto const omega = std::sqrt(this->c2 / me); // n = 1, |q| = 1

    this->species.push_back({+1., 1e6, 0., typename TestFixture::ParticleArray_t{box}});
    this->species.push_back({-1., me, 0., typename TestFixture::ParticleArray_t{box}});
    for (int cell = first; cell <= last; ++cell)
        for (int i = 0; i < ppc; ++i)
        {
            double const delta = (i + 0.5) / ppc;
            double const x     = (cell + delta) * dx;
            double const vx    = v0 * std::sin(std::numbers::pi * (x - first * dx) / slab);
            this->species[0].particles.emplace_back(
                Particle<1>{1. / ppc, {cell}, {delta}, {0., 0., 0.}});
            this->species[1].particles.emplace_back(
                Particle<1>{1. / ppc, {cell}, {delta}, {vx, 0., 0.}});
        }

    // the electron in the middle of the slab
    auto const probe = static_cast<std::size_t>((last + 1 - first) * ppc / 2);
    ASSERT_GT(this->species[1].particles[probe].v[0], 0.9 * v0);

    this->startUp();

    std::vector<double> crossings;
    double previous = this->species[1].particles[probe].v[0];
    for (std::size_t n = 1; n <= 120; ++n)
    {
        this->step();
        auto const v = this->species[1].particles[probe].v[0];
        if (v * previous < 0)
            crossings.push_back(dt * (n - 1 + previous / (previous - v)));
        previous = v;
    }
    ASSERT_GE(crossings.size(), 6u);

    auto const halfPeriod  = (crossings.back() - crossings.front()) / (crossings.size() - 1);
    auto const measured    = std::numbers::pi / halfPeriod;
    auto const leapfrogged = 2 / dt * std::asin(omega * dt / 2);
    std::cout << "plasma oscillation: measured " << measured << " expected " << leapfrogged << "\n";
    EXPECT_NEAR(measured / leapfrogged, 1., 5e-3);
}


// no particles: a standing light wave B_z = cos(kx) cos(wt) follows the Yee dispersion
// sin(w dt / 2) / (c dt) = sin(k dx / 2) / dx. Checked in the middle of the box, before the
// stale ghost nodes can be felt there.
TYPED_TEST(FKLeapfrog1D, lightWaveFollowsTheYeeDispersion)
{
    auto const dx = this->dx;
    auto const dt = this->dt;
    auto const c  = std::sqrt(this->c2);
    auto const k  = 2 * std::numbers::pi / (this->nCells * dx) * 4;
    auto& Bz      = this->B(Component::Z);
    this->fill_(Bz, [&](double x) { return std::cos(k * x); });

    this->startUp();

    auto const omega = 2 / dt * std::asin(c * dt / dx * std::sin(k * dx / 2));
    auto const mid   = this->layout.physicalStartIndex(QtyCentering::dual, Direction::X)
                     + this->nCells / 2; // Bz is dual in 1D
    auto const xmid = (this->nCells / 2 + 0.5) * dx;

    for (std::size_t n = 1; n <= 20; ++n)
    {
        this->step();
        auto const expected = std::cos(k * xmid) * std::cos(omega * n * dt);
        EXPECT_NEAR(Bz(mid), expected, 2e-3) << "at step " << n;
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
