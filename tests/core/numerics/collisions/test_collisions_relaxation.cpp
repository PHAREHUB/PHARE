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
double isotropizationRatio(double const dt, bool const randomWeights)
{
    auto const nu0 = 0.1, mass = 1.;
    std::mt19937_64 rng{5};
    TestIons ions;
    ions.populations.emplace_back(mass, 1., nbrCells);
    ions.populations[0].load(2000, 1., {0., 0., 0.}, {0.2, 0.3, 0.3}, rng, randomWeights);
    Collisions<TestIons> collide{{nu0, 1, ScatteringLaw::nanbu, seed}};

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
    EXPECT_NEAR(isotropizationRatio(0.0025, false), 1., 0.1);
}

TEST(CollisionRelaxation, temperatureAnisotropyRelaxesAtTheNRLRateWithUnequalWeights)
{
    EXPECT_NEAR(isotropizationRatio(0.0025, true), 1., 0.1);
}

TEST(CollisionRelaxation, unequalWeightsRateConvergesInDt)
{
    auto const coarse = isotropizationRatio(0.01, true);
    auto const fine   = isotropizationRatio(0.0025, true);
    EXPECT_LT(std::abs(1. - fine), std::abs(1. - coarse));
}


TEST(CollisionRelaxation, protonsAndAlphasEquilibrateAtTheNRLRate)
{
    auto const nu0 = 0.01, dt = 0.0025;
    std::mt19937_64 rng{11};
    TestIons ions;
    ions.populations.emplace_back(1., 1., nbrCells);
    ions.populations.emplace_back(4., 2., nbrCells);
    ions.populations[0].load(5000, 1., {0., 0., 0.}, {0.2, 0.2, 0.2}, rng, false);
    ions.populations[1].load(2500, 0.25, {0., 0., 0.}, {0.05, 0.05, 0.05}, rng, false);
    Collisions<TestIons> collide{{nu0, 1, ScatteringLaw::nanbu, seed}};

    auto temperatures = [&]() {
        return std::array{Moments{ions.populations[0], nbrCells}.isotropic(),
                          Moments{ions.populations[1], nbrCells}.isotropic()};
    };
    auto const T0  = temperatures();
    auto T         = T0;
    auto const nH  = Moments{ions.populations[0], nbrCells}.density;
    auto const nHe = Moments{ions.populations[1], nbrCells}.density;

    auto const rate = [&](std::array<double, 2> const& t) {
        return std::array{equilibrationRate(nu0, 1., 2., nHe, 1., 4., t[0], t[1]) * (t[1] - t[0]),
                          equilibrationRate(nu0, 2., 1., nH, 4., 1., t[1], t[0]) * (t[0] - t[1])};
    };

    // about a tenth of the temperature difference relaxes
    auto constexpr steps = 52;
    int step             = 0;
    advance(collide, ions, steps, step, dt);
    for (int i = 0; i < 10 * steps; ++i)
        T = rk4(T, dt / 10., rate);

    auto const measured  = temperatures();
    auto const reference = std::log((T[0] - T[1]) / (T0[0] - T0[1]));
    EXPECT_LT(reference, -0.05);
    EXPECT_NEAR(std::log((measured[0] - measured[1]) / (T0[0] - T0[1])) / reference, 1., 0.1);
}


TEST(CollisionRelaxation, equalWeightsConserveMomentumAndEnergyExactly)
{
    std::mt19937_64 rng{17};
    TestIons ions;
    ions.populations.emplace_back(1., 1., nbrCells);
    ions.populations.emplace_back(4., 2., nbrCells);
    ions.populations[0].load(200, 1., {0.1, 0., 0.}, {0.2, 0.3, 0.3}, rng, false);
    ions.populations[1].load(200, 1., {-0.1, 0., 0.}, {0.1, 0.1, 0.1}, rng, false);
    Collisions<TestIons> collide{{0.05, 1, ScatteringLaw::nanbu, seed}};

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
    Collisions<TestIons> collide{{0.1, 1, ScatteringLaw::nanbu, seed}};

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
