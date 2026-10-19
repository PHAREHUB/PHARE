#include "core/numerics/collisions/collisions.hpp"

#include "test_collisions_fixtures.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <set>

using namespace PHARE::core;
using namespace PHARE::core::collisions;



TEST(CollisionRandom, splitmix64MatchesTheReferenceGeneratorSeededWithZero)
{
    // the reference SplitMix64 seeded with 0 outputs splitmix64(0), splitmix64(golden), ...
    EXPECT_EQ(splitmix64(0), 0xe220a8397b1dcdafULL);
    EXPECT_EQ(splitmix64(0x9e3779b97f4a7c15ULL), 0x6e789e6aa1b965f4ULL);
}

TEST(CollisionRandom, unitDrawsAreInZeroOneClosedAtOne)
{
    EXPECT_GT(toUnit(0), 0.);
    EXPECT_EQ(toUnit(~std::uint64_t{0}), 1.);
}

TEST(CollisionRandom, neighbourCellsHaveUnrelatedStreams)
{
    // a stream must not be a shifted copy of its neighbour's
    auto const key = std::uint64_t{12345};
    std::set<double> draws;
    for (std::uint64_t cell = 0; cell < 64; ++cell)
    {
        auto stream = Stream{combine(key, cell)};
        for (int i = 0; i < 64; ++i)
            draws.insert(stream());
    }
    EXPECT_EQ(draws.size(), 64u * 64u);
}



TEST(CollisionKinematics, scatteringKeepsPairMomentumEnergyAndRelativeSpeed)
{
    std::mt19937_64 rng{1};
    std::uniform_real_distribution<double> uniform{-1., 1.};

    for (int trial = 0; trial < 1000; ++trial)
    {
        std::array<double, 3> va{uniform(rng), uniform(rng), uniform(rng)};
        std::array<double, 3> vb{uniform(rng), uniform(rng), uniform(rng)};
        auto const ma     = 1. + 3. * (uniform(rng) + 1.);
        auto const mb     = 1. + 3. * (uniform(rng) + 1.);
        auto const cosChi = uniform(rng);
        auto const phi    = M_PI * (uniform(rng) + 1.);

        auto const p0 = [&](std::size_t i) { return ma * va[i] + mb * vb[i]; };
        auto const e  = [&]() {
            double sum = 0;
            for (std::size_t i = 0; i < 3; ++i)
                sum += 0.5 * (ma * va[i] * va[i] + mb * vb[i] * vb[i]);
            return sum;
        };
        auto const u = [&]() { return std::array{va[0] - vb[0], va[1] - vb[1], va[2] - vb[2]}; };

        std::array const p{p0(0), p0(1), p0(2)};
        auto const e0 = e();
        auto const u0 = u();

        scatter(va, vb, ma, mb, cosChi, phi, true, true);

        auto const u1 = u();
        for (std::size_t i = 0; i < 3; ++i)
            EXPECT_NEAR(p0(i), p[i], 1e-13);
        EXPECT_NEAR(e(), e0, 1e-13);
        auto const dot  = u0[0] * u1[0] + u0[1] * u1[1] + u0[2] * u1[2];
        auto const norm = [](auto const& x) { return std::hypot(x[0], x[1], x[2]); };
        EXPECT_NEAR(norm(u1), norm(u0), 1e-13);
        EXPECT_NEAR(dot / (norm(u0) * norm(u1)), cosChi, 1e-12);
    }
}

TEST(CollisionKinematics, relativeVelocityAlongZIsDeflectedByChi)
{
    std::array<double, 3> va{0., 0., 1.}, vb{0., 0., -1.};
    scatter(va, vb, 1., 1., 0., 0.3, true, true);

    // u was (0, 0, 2), chi = pi/2 turns it into the x-y plane, momentum stays zero
    EXPECT_NEAR(va[2] - vb[2], 0., 1e-14);
    EXPECT_NEAR(std::hypot(va[0] - vb[0], va[1] - vb[1]), 2., 1e-14);
    for (std::size_t i = 0; i < 3; ++i)
        EXPECT_NEAR(va[i] + vb[i], 0., 1e-14);
}



// mean of cos chi for the law as implemented, against which samples are tested
double lawMeanCosChi(double const s)
{
    if (s < 0.1)
        return 1. - s; // cos chi = 1 + s ln U, truncated at -1 with negligible weight
    auto const A
        = s < 3. ? 1.
                       / (0.0056958 + 0.9560202 * s - 0.508139 * s * s + 0.47913906 * s * s * s
                          - 0.12788975 * s * s * s * s + 0.02389567 * s * s * s * s * s)
                 : 3. * std::exp(-s);
    return s < 6. ? 1. / std::tanh(A) - 1. / A : 0.;
}

TEST(CollisionAngle, nanbuSamplesTheFittedLawWhichFollowsExpMinusS)
{
    for (auto const s : {0.01, 0.05, 0.3, 1., 2.5, 4., 8.})
    {
        auto draw        = Stream{combine(7, std::bit_cast<std::uint64_t>(s))};
        auto constexpr n = 200000;
        double sum = 0, sq = 0;
        for (int i = 0; i < n; ++i)
        {
            auto const c = cosChi<ScatteringLaw::nanbu>(s, draw);
            ASSERT_GE(c, -1.);
            ASSERT_LE(c, 1.);
            sum += c;
            sq += c * c;
        }
        auto const mean  = sum / n;
        auto const sigma = std::sqrt((sq / n - mean * mean) / n);
        EXPECT_NEAR(mean, lawMeanCosChi(s), 5. * sigma + 1e-12) << "s = " << s;

        // Nanbu's law has <cos chi> = exp(-s); the fit, the small-s form and the isotropic
        // cut-off beyond s = 6 approximate it
        auto const tolerance = s < 0.1 ? 0.5 * s * s : s < 6. ? 0.01 * std::exp(-s) : std::exp(-6.);
        EXPECT_NEAR(lawMeanCosChi(s), std::exp(-s), tolerance) << "s = " << s;
    }
}

TEST(CollisionAngle, takizukaAbeMeanDeflectionIsSAtSmallS)
{
    auto const s     = 0.01;
    auto draw        = Stream{99};
    auto constexpr n = 400000;
    double sum       = 0;
    for (int i = 0; i < n; ++i)
        sum += 1. - cosChi<ScatteringLaw::takizuka_abe>(s, draw);
    EXPECT_NEAR(sum / n / s, 1., 0.02);
}

TEST(CollisionAngle, bobylevIsDeterministic)
{
    auto draw = Stream{0};
    EXPECT_DOUBLE_EQ(cosChi<ScatteringLaw::bobylev>(0.3, draw), 0.7);
    EXPECT_DOUBLE_EQ(cosChi<ScatteringLaw::bobylev>(5., draw), -1.);
    EXPECT_EQ(draw.index, 0u);
}



TEST(CollisionPairing, intraPairsUseEveryParticleOnceOrTheOddTriangleTwiceAtHalfStrength)
{
    for (std::size_t n = 2; n < 12; ++n)
    {
        std::map<std::size_t, double> strength;
        forEachIntraPair(n, [&](auto const i, auto const j, double const factor) {
            EXPECT_NE(i, j);
            strength[i] += factor;
            strength[j] += factor;
        });
        ASSERT_EQ(strength.size(), n);
        for (auto const& [i, total] : strength)
            EXPECT_DOUBLE_EQ(total, 1.) << "n = " << n << " i = " << i;
    }
}

TEST(CollisionPairing, interPairsMeetEveryLargeOnceAndShareTheSmallEvenly)
{
    for (std::size_t nLarge = 1; nLarge < 20; ++nLarge)
        for (std::size_t nSmall = 1; nSmall <= nLarge; ++nSmall)
        {
            std::map<std::size_t, std::size_t> large, small;
            forEachInterPair(nLarge, nSmall, [&](auto const i, auto const j) {
                ++large[i];
                ++small[j];
            });
            EXPECT_EQ(large.size(), nLarge);
            for (auto const& [_, count] : large)
                EXPECT_EQ(count, 1u);
            EXPECT_EQ(small.size(), nSmall);
            for (auto const& [_, count] : small)
            {
                EXPECT_GE(count, nLarge / nSmall);
                EXPECT_LE(count, nLarge / nSmall + 1);
            }
        }
}



// Each pair adds strength x larger weight x acceptance to the scattering of each partner.
// Averaged over shuffles, a particle must see the density of the particles it can scatter on:
// the others of its group, skipping its siblings (no relative velocity), or the other group.
TEST(CollisionPairing, intraPairsGiveEachParticleTheDensityOfTheOthers)
{
    std::mt19937_64 rng{9};
    std::uniform_real_distribution<double> uniform{0.01, 1.};
    auto constexpr shuffles = 20000;

    for (std::size_t const n : {2u, 3u, 4u, 5u, 7u})
    {
        std::vector<double> weight(n);
        std::vector<std::size_t> sibling(n), order(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            weight[i]  = uniform(rng);
            sibling[i] = order[i] = i;
        }
        if (n >= 4) // children of one split parent: same weight, same velocity
        {
            weight[1]  = weight[0];
            sibling[1] = 0;
        }

        std::vector<double> seen(n, 0.);
        for (int shuffle = 0; shuffle < shuffles; ++shuffle)
        {
            std::shuffle(order.begin(), order.end(), rng);
            forEachIntraPair(n, [&](auto const i, auto const j, double const factor) {
                auto const a = order[i], b = order[j];
                if (sibling[a] == sibling[b])
                    return;
                auto const [pA, pB] = acceptance(weight[a], weight[b]);
                auto const strength = intraStrength(n, factor) * std::max(weight[a], weight[b]);
                seen[a] += strength * pA;
                seen[b] += strength * pB;
            });
        }

        for (std::size_t i = 0; i < n; ++i)
        {
            double expected = 0;
            for (std::size_t j = 0; j < n; ++j)
                if (sibling[j] != sibling[i])
                    expected += weight[j];
            EXPECT_NEAR(seen[i] / shuffles, expected, 0.03 * expected) << "n " << n << " i " << i;
        }
    }
}

TEST(CollisionPairing, interPairsGiveEachParticleTheDensityOfTheOtherGroup)
{
    std::mt19937_64 rng{10};
    std::uniform_real_distribution<double> uniform{0.01, 1.};
    auto constexpr shuffles = 20000;

    for (auto const [nLarge, nSmall] : {std::pair{5u, 2u}, {7u, 3u}, {4u, 4u}, {6u, 1u}})
    {
        std::vector<double> wLarge(nLarge), wSmall(nSmall);
        std::vector<std::size_t> oLarge(nLarge), oSmall(nSmall);
        for (std::size_t i = 0; i < nLarge; ++i)
            wLarge[i] = uniform(rng), oLarge[i] = i;
        for (std::size_t i = 0; i < nSmall; ++i)
            wSmall[i] = uniform(rng), oSmall[i] = i;

        std::vector<double> seenLarge(nLarge, 0.), seenSmall(nSmall, 0.);
        for (int shuffle = 0; shuffle < shuffles; ++shuffle)
        {
            std::shuffle(oLarge.begin(), oLarge.end(), rng);
            std::shuffle(oSmall.begin(), oSmall.end(), rng);
            forEachInterPair(nLarge, nSmall, [&](auto const i, auto const j) {
                auto const a = oLarge[i], b = oSmall[j];
                auto const [pA, pB] = acceptance(wLarge[a], wSmall[b]);
                auto const strength = interStrength(nSmall) * std::max(wLarge[a], wSmall[b]);
                seenLarge[a] += strength * pA;
                seenSmall[b] += strength * pB;
            });
        }

        auto const nL = std::accumulate(wLarge.begin(), wLarge.end(), 0.);
        auto const nS = std::accumulate(wSmall.begin(), wSmall.end(), 0.);
        for (auto const seen : seenLarge)
            EXPECT_NEAR(seen / shuffles, nS, 0.03 * nS);
        for (auto const seen : seenSmall)
            EXPECT_NEAR(seen / shuffles, nL, 0.03 * nL);
    }
}



TEST(CollisionOperator, isReproducibleAndIndependentOfStorageOrder)
{
    auto const nbrCells = 4;

    auto makeIons = [&](bool reversed) {
        TestIons ions;
        ions.populations.emplace_back(1., 1., nbrCells);
        std::mt19937_64 local{3};
        auto& pop = ions.populations.back();
        TestPopulation tmp{1., 1., nbrCells};
        tmp.load(50, 1., {0., 0., 0.}, {0.2, 0.3, 0.3}, local, true);
        auto particles = tmp.domainParticles().vector();
        if (reversed)
            std::reverse(particles.begin(), particles.end());
        for (auto const& p : particles)
            pop.domainParticles().push_back(p);
        return ions;
    };

    auto a = makeIons(false);
    auto b = makeIons(true);
    Collisions<TestIons> collideA{{1e-2, 1, ScatteringLaw::nanbu, 42}};
    Collisions<TestIons> collideB{{1e-2, 1, ScatteringLaw::nanbu, 42}};
    collideA(a, Box<int, 1>{{0}, {nbrCells - 1}}, 1234, 0.01);
    collideB(b, Box<int, 1>{{0}, {nbrCells - 1}}, 1234, 0.01);

    auto sorted = [](auto const& pop) {
        auto v = pop.domainParticles().vector();
        std::sort(v.begin(), v.end(), [](auto const& x, auto const& y) {
            return std::tie(x.iCell, x.delta, x.weight) < std::tie(y.iCell, y.delta, y.weight);
        });
        return v;
    };
    auto const pa = sorted(a.populations[0]);
    auto const pb = sorted(b.populations[0]);
    ASSERT_EQ(pa.size(), pb.size());
    for (std::size_t i = 0; i < pa.size(); ++i)
        EXPECT_EQ(pa[i].v, pb[i].v);
}

TEST(CollisionOperator, siblingsWithEqualVelocitiesAreSkippedWithoutNaN)
{
    TestIons ions;
    ions.populations.emplace_back(1., 1., 1);
    Particle<1> particle;
    particle.weight = 0.1;
    particle.iCell  = {0};
    particle.v      = {0.1, -0.2, 0.3};
    for (int i = 0; i < 4; ++i)
    {
        particle.delta = {0.1 + 0.2 * i};
        ions.populations[0].domainParticles().push_back(particle);
    }

    Collisions<TestIons> collide{{1., 1, ScatteringLaw::nanbu, 1}};
    collide(ions, Box<int, 1>{{0}, {0}}, 7, 0.1);

    EXPECT_EQ(collide.counters().pairs, 2u);
    EXPECT_EQ(collide.counters().skipped, 2u);
    for (auto const& p : ions.populations[0].domainParticles())
        EXPECT_EQ(p.v, particle.v);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
