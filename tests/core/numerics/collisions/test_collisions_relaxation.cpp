// relaxation of homogeneous plasmas under collisions only, against the NRL formulary rates,
// in PHARE units (nu0 = Coulomb logarithm times e^4 n0 / (4 pi eps0^2 mp^2 VA^3 Omega_p))

#include "core/numerics/collisions/collisions.hpp"

#include "test_collisions_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <numbers>
#include <random>

using namespace PHARE::core;

namespace
{
auto constexpr nbrCells = 20;
auto const cells        = Box<int, 1>{{0}, {nbrCells - 1}};
auto constexpr seed     = std::uint64_t{2026};


void advance(Collisions<TestIons>& collide, TestIons& ions, int const steps, int& step,
             double const dt)
{
    for (auto const end = step + steps; step < end; ++step)
        collide(ions, cells, collisions::combine(seed, static_cast<std::uint64_t>(step)), dt);
}


// NRL: dT_perp/dt = -1/2 dT_par/dt = -nu_T (T_perp - T_par)
struct Isotropization
{
    double nu0, density, mass;

    std::array<double, 2> operator()(std::array<double, 2> const& T) const // {T_par, T_perp}
    {
        auto const A = T[1] / T[0] - 1.;
        auto const F = [&]() {
            if (std::abs(A) < 1e-6)
                return 4. / 15.;
            auto const r = std::sqrt(std::abs(A));
            auto const g = A > 0 ? std::atan(r) / r : std::atanh(r) / r;
            return (-3. + (A + 3.) * g) / (A * A);
        }();
        auto const nuT
            = nu0 * density
              / (2. * std::sqrt(std::numbers::pi) * std::sqrt(mass) * std::pow(T[0], 1.5)) * F;
        auto const dTperp = -nuT * (T[1] - T[0]);
        return {-2. * dTperp, dTperp};
    }
};


// ln of the relaxed fraction of the anisotropy, simulated over reference, at an early time
// (NRL assumes bi-maxwellians, true only at first). Binary collisions converge to the NRL rate
// as dt -> 0, unequal weights roughly doubling the dt error (Wang et al. 2008 for the first).
double isotropizationRatio(double const dt, bool const randomWeights, std::uint64_t const rngSeed)
{
    auto const nu0 = 0.1, mass = 1.;
    std::mt19937_64 rng{rngSeed};
    TestIons ions;
    ions.populations.emplace_back(mass, 1., nbrCells);
    ions.populations[0].load(2000, 1., {0., 0., 0.}, {0.2, 0.3, 0.3}, rng, randomWeights);
    Collisions<TestIons> collide{{nu0, 1, ScatteringLaw::nanbu, false, seed}};

    auto const m0 = Moments{ions.populations[0], nbrCells};
    std::array<double, 2> T{m0.parallel(), m0.perpendicular()};
    auto const dT0  = T[1] - T[0];
    auto const rate = Isotropization{nu0, m0.density, mass};

    // about a quarter of the anisotropy relaxes
    auto const initialRate = -3. * rate(T)[1] / dT0;
    auto const steps       = static_cast<int>(0.3 / initialRate / dt);
    int step               = 0;
    advance(collide, ions, steps, step, dt);
    for (int i = 0; i < 10 * steps; ++i)
        T = rk4(T, dt / 10., rate);

    auto const m = Moments{ions.populations[0], nbrCells};
    return std::log((m.perpendicular() - m.parallel()) / dT0) / std::log((T[1] - T[0]) / dT0);
}


// the ratio fluctuates by about 0.04 from one initial sampling to another
double meanIsotropizationRatio(double const dt, bool const randomWeights)
{
    auto constexpr samplings = 4;
    double sum               = 0;
    for (std::uint64_t rngSeed = 1; rngSeed <= samplings; ++rngSeed)
        sum += isotropizationRatio(dt, randomWeights, rngSeed);
    return sum / samplings;
}


// NRL: dT_a/dt = nu_ab (T_b - T_a), nu_ab = 2/3 sqrt(2/pi) nu0 Za^2 Zb^2 n_b sqrt(ma mb)
//                                           / (ma Tb + mb Ta)^(3/2)
double equilibrationRate(double const nu0, double const za, double const zb, double const nb,
                         double const ma, double const mb, double const Ta, double const Tb)
{
    return 2. / 3. * std::sqrt(2. / std::numbers::pi) * nu0 * za * za * zb * zb * nb
           * std::sqrt(ma * mb) / std::pow(ma * Tb + mb * Ta, 1.5);
}

} // namespace



TEST(CollisionRelaxation, temperatureAnisotropyRelaxesAtTheNRLRate)
{
    EXPECT_NEAR(meanIsotropizationRatio(0.0025, false), 1., 0.1);
}

TEST(CollisionRelaxation, temperatureAnisotropyRelaxesAtTheNRLRateWithUnequalWeights)
{
    // weights over two decades double the dt error, measured 0.947 +- 0.019 here
    EXPECT_NEAR(meanIsotropizationRatio(0.0025, true), 1., 0.11);
}

TEST(CollisionRelaxation, unequalWeightsRateConvergesInDt)
{
    auto const coarse = meanIsotropizationRatio(0.01, true);
    auto const fine   = meanIsotropizationRatio(0.0025, true);
    EXPECT_LT(std::abs(1. - fine), std::abs(1. - coarse));
}


namespace
{
// protons (n = 1) and alphas (n = 0.25), between groups only since each is a maxwellian
struct ProtonsAndAlphas
{
    double nu0;
    TestIons ions;
    Collisions<TestIons> collide;
    int step = 0;

    ProtonsAndAlphas(double const nu0_, std::array<double, 3> const& protonDrift,
                     double const protonVth, std::array<double, 3> const& alphaDrift,
                     double const alphaVth, bool const randomWeights, bool const exact)
        : nu0{nu0_}
        , collide{{nu0_, 1, ScatteringLaw::nanbu, exact, seed}}
    {
        std::mt19937_64 rng{11};
        ions.populations.emplace_back(1., 1., nbrCells);
        ions.populations.emplace_back(4., 2., nbrCells);
        auto const iso = [](double v) { return std::array{v, v, v}; };
        ions.populations[0].load(5000, 1., protonDrift, iso(protonVth), rng, randomWeights);
        ions.populations[1].load(2500, 0.25, alphaDrift, iso(alphaVth), rng, randomWeights);
    }

    auto moments() const
    {
        return std::array{Moments{ions.populations[0], nbrCells},
                          Moments{ions.populations[1], nbrCells}};
    }

    void advance(int const steps, double const dt) { ::advance(collide, ions, steps, step, dt); }
};


// ln of the relaxed fraction of the temperature difference, simulated over NRL, early time
double equilibrationRatio(bool const randomWeights, bool const exact = false)
{
    auto const nu0 = 0.01, dt = 0.0025;
    ProtonsAndAlphas plasma{nu0, {0., 0., 0.}, 0.2, {0., 0., 0.}, 0.05, randomWeights, exact};

    auto const m0 = plasma.moments();
    std::array const T0{m0[0].isotropic(), m0[1].isotropic()};
    auto T = T0;

    auto const rate = [&](std::array<double, 2> const& t) {
        return std::array{
            equilibrationRate(nu0, 1., 2., m0[1].density, 1., 4., t[0], t[1]) * (t[1] - t[0]),
            equilibrationRate(nu0, 2., 1., m0[0].density, 4., 1., t[1], t[0]) * (t[0] - t[1])};
    };

    // about a tenth of the temperature difference relaxes
    auto constexpr steps = 52;
    plasma.advance(steps, dt);
    for (int i = 0; i < 10 * steps; ++i)
        T = rk4(T, dt / 10., rate);

    auto const m         = plasma.moments();
    auto const reference = std::log((T[0] - T[1]) / (T0[0] - T0[1]));
    EXPECT_LT(reference, -0.05);
    return std::log((m[0].isotropic() - m[1].isotropic()) / (T0[0] - T0[1])) / reference;
}


// friction between maxwellians at the same temperature with a small relative drift:
// dU_a/dt = nu_ab (U_b - U_a), the momentum transfer frequency being
// nu_ab = nu_eps^{a|b} (ma + mb) / (2 ma) (Burgers 1969, Schunk 1977)
double frictionRatio(bool const randomWeights, bool const exact = false, double const boost = 0.)
{
    auto const nu0 = 0.01, dt = 0.0025;
    // T = 0.04 for both, zero total momentum, drift small against the relative thermal speed;
    // `boost` moves the whole plasma, which must not change anything
    ProtonsAndAlphas plasma{nu0, {boost + 0.03, 0., 0.}, 0.2,  {boost - 0.03, 0., 0.},
                            0.1, randomWeights,          exact};

    auto const m0 = plasma.moments();
    auto const T  = m0[0].isotropic();
    auto const nuPA
        = equilibrationRate(nu0, 1., 2., m0[1].density, 1., 4., T, T) * (1. + 4.) / (2. * 1.);
    auto const nuAP
        = equilibrationRate(nu0, 2., 1., m0[0].density, 4., 1., T, T) * (4. + 1.) / (2. * 4.);
    auto const dU0 = m0[0].velocity[0] - m0[1].velocity[0];

    // about a fifth of the relative drift is lost
    auto const steps = static_cast<int>(0.2 / (nuPA + nuAP) / dt);
    plasma.advance(steps, dt);

    auto const m = plasma.moments();
    return std::log((m[0].velocity[0] - m[1].velocity[0]) / dU0) / (-(nuPA + nuAP) * steps * dt);
}
} // namespace


TEST(CollisionRelaxation, protonsAndAlphasEquilibrateAtTheNRLRate)
{
    EXPECT_NEAR(equilibrationRatio(false), 1., 0.1);
}

TEST(CollisionRelaxation, protonsAndAlphasEquilibrateAtTheNRLRateWithUnequalWeights)
{
    EXPECT_NEAR(equilibrationRatio(true), 1., 0.1);
}

TEST(CollisionRelaxation, protonsAndAlphasDriftRelaxesAtTheMomentumTransferRate)
{
    EXPECT_NEAR(frictionRatio(false), 1., 0.15);
}

TEST(CollisionRelaxation, protonsAndAlphasDriftRelaxesAtTheMomentumTransferRateWithUnequalWeights)
{
    EXPECT_NEAR(frictionRatio(true), 1., 0.15);
}


// with the momentum and energy fix, which must keep the exchanges between groups
TEST(CollisionRelaxation, protonsAndAlphasEquilibrateAtTheNRLRateWithTheConservationFix)
{
    EXPECT_NEAR(equilibrationRatio(true, true), 1., 0.1);
}

TEST(CollisionRelaxation, protonsAndAlphasDriftRelaxesAtTheMomentumTransferRateWithTheFix)
{
    EXPECT_NEAR(frictionRatio(true, true), 1., 0.15);
}

TEST(CollisionRelaxation, aBoostedPlasmaRelaxesTheSameWithTheFix)
{
    EXPECT_NEAR(frictionRatio(true, true, 2.), 1., 0.15);
}


// Angus et al. 2025, Test 1: a beam of 1/10 of the density, macro-particles ten times lighter,
// in one species. It must end at rest in the frame of the centre of mass and at the temperature
// given by the conservation of energy
TEST(CollisionRelaxation, beamAndBackgroundReachTheConservedEndState)
{
    auto constexpr cells = 4;
    std::mt19937_64 rng{31};
    TestIons ions;
    ions.populations.emplace_back(1., 1., cells);
    ions.populations.emplace_back(1., 1., cells); // the beam, same species
    ions.populations[0].load(1000, 1., {0., 0., 0.}, {0.1, 0.1, 0.1}, rng, false);
    ions.populations[1].load(1000, 0.1, {0.5, 0., 0.}, {0.1, 0.1, 0.1}, rng, false);
    Collisions<TestIons> collide{{1., 1, ScatteringLaw::nanbu, true, seed}};

    auto const m0
        = std::array{Moments{ions.populations[0], cells}, Moments{ions.populations[1], cells}};
    auto const n  = m0[0].density + m0[1].density;
    auto const Uf = (m0[0].density * m0[0].velocity[0] + m0[1].density * m0[1].velocity[0]) / n;
    auto const kinetic = [&](auto const& m) {
        return 0.5 * m.density
               * (std::pow(m.velocity[0], 2) + std::pow(m.velocity[1], 2)
                  + std::pow(m.velocity[2], 2));
    };
    auto const Tf = (1.5 * (m0[0].density * m0[0].isotropic() + m0[1].density * m0[1].isotropic())
                     + kinetic(m0[0]) + kinetic(m0[1]) - 0.5 * n * Uf * Uf)
                    / (1.5 * n);

    for (int step = 0; step < 400; ++step)
        collide(ions, Box<int, 1>{{0}, {cells - 1}}, collisions::combine(seed, step), 0.01);
    EXPECT_EQ(collide.counters().unconverged, 0u);

    for (auto const& pop : ions)
    {
        auto const m = Moments{pop, cells};
        EXPECT_NEAR(m.velocity[0], Uf, 0.03 * Uf);
        EXPECT_NEAR(m.isotropic(), Tf, 0.05 * Tf);
    }
}


TEST(CollisionRelaxation, equalWeightsConserveMomentumAndEnergyExactly)
{
    std::mt19937_64 rng{17};
    TestIons ions;
    ions.populations.emplace_back(1., 1., nbrCells);
    ions.populations.emplace_back(4., 2., nbrCells);
    ions.populations[0].load(200, 1., {0.1, 0., 0.}, {0.2, 0.3, 0.3}, rng, false);
    ions.populations[1].load(200, 1., {-0.1, 0., 0.}, {0.1, 0.1, 0.1}, rng, false);
    Collisions<TestIons> collide{{0.05, 1, ScatteringLaw::nanbu, false, seed}};

    auto totals = [&]() {
        std::array<double, 4> sum{};
        for (auto const& pop : ions)
            for (auto const& p : pop.domainParticles())
            {
                for (std::size_t i = 0; i < 3; ++i)
                    sum[i] += pop.mass() * p.weight * p.v[i];
                sum[3] += 0.5 * pop.mass() * p.weight
                          * (p.v[0] * p.v[0] + p.v[1] * p.v[1] + p.v[2] * p.v[2]);
            }
        return sum;
    };

    auto const before = totals();
    int step          = 0;
    advance(collide, ions, 50, step, 0.01);
    auto const after = totals();
    for (std::size_t i = 0; i < 4; ++i)
        EXPECT_NEAR(after[i], before[i], 1e-12 * std::abs(before[3]));
}


TEST(CollisionRelaxation, maxwellianStaysMaxwellianWithUnequalWeights)
{
    std::mt19937_64 rng{23};
    TestIons ions;
    ions.populations.emplace_back(1., 1., nbrCells);
    ions.populations[0].load(1000, 1., {0., 0., 0.}, {0.2, 0.2, 0.2}, rng, true);
    Collisions<TestIons> collide{{0.1, 1, ScatteringLaw::nanbu, false, seed}};

    auto const T0 = Moments{ions.populations[0], nbrCells}.isotropic();
    int step      = 0;
    advance(collide, ions, 200, step, 0.01);
    auto const m = Moments{ions.populations[0], nbrCells};

    // weighted scattering conserves only on average: drift of a few percent of the noise
    EXPECT_NEAR(m.isotropic(), T0, 0.02 * T0);
    EXPECT_NEAR(m.parallel(), m.perpendicular(), 0.03 * T0);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
