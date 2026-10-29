
#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/data/ions/particle_initializers/maxwellian_particle_initializer.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_utilities.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/random.hpp"

#include "gtest/gtest.h"

#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <optional>

#include "tests/initializer/init_functions.hpp"
using namespace PHARE::initializer::test_fn::func_1d; // density/etc are here

using namespace PHARE::core;
using namespace PHARE::initializer;


class AMaxwellianParticleInitializer1D : public ::testing::Test
{
private:
    using GridLayoutT       = PHARE::core::PHARE_Types<PHARE::SimOpts{1, 1}>::Hybrid::GridLayout_t;
    using ParticleArrayT    = ParticleArray<1>;
    using InitFunctionArray = std::array<InitFunction<1>, 3>;

public:
    AMaxwellianParticleInitializer1D()
        : layout{{{0.1}}, {{50}}, Point{0.}, Box{Point{50}, Point{99}}}
        , particles{layout.AMRBox()}
        , initializer{std::make_unique<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>>(
              density, InitFunctionArray{vx, vy, vz}, InitFunctionArray{vthx, vthy, vthz},
              nbrParticlesPerCell)}
    {
        //
    }


    GridLayoutT layout;
    ParticleArrayT particles;
    std::uint32_t nbrParticlesPerCell{10000};
    std::unique_ptr<MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT>> initializer;
};




TEST_F(AMaxwellianParticleInitializer1D, loadsTheCorrectNbrOfParticles)
{
    auto nbrCells             = layout.nbrCells();
    auto expectedNbrParticles = nbrParticlesPerCell * nbrCells[0];
    initializer->loadParticles(particles, layout);
    EXPECT_EQ(expectedNbrParticles, particles.size());
}




TEST_F(AMaxwellianParticleInitializer1D, loadsParticlesInTheDomain)
{
    initializer->loadParticles(particles, layout);
    for (auto const& particle : particles)
    {
        EXPECT_TRUE(particle.iCell[0] >= 50 && particle.iCell[0] <= 99);
        auto pos       = positionAsPoint(particle, layout);
        auto endDomain = (layout.AMRBox().upper[0] + 1) * layout.meshSize()[0];

        if (!((pos[0] > 0.) and (pos[0] < endDomain)))
            std::cout << "position : " << pos[0] << " not in domain (0," << endDomain << ")\n";
        EXPECT_TRUE(pos[0] > 0. && pos[0] < endDomain);
    }
}




// per-cell seeding: particles of a cell depend on (seed, population, level, AMR cell) only

template<std::size_t dim>
struct SeededLoader
{
    using GridLayoutT =
        typename PHARE::core::PHARE_Types<PHARE::SimOpts{dim, 1}>::Hybrid::GridLayout_t;
    using ParticleArrayT    = ParticleArray<dim>;
    using InitFunctionArray = std::array<InitFunction<dim>, 3>;
    using Cell              = std::array<int, dim>;
    using CellParticles     = std::map<Cell, std::vector<Particle<dim>>>;

    static auto functions()
    {
        if constexpr (dim == 1)
            return std::make_tuple(InitFunction<1>{density}, InitFunctionArray{vx, vy, vz},
                                   InitFunctionArray{vthx, vthy, vthz});
        else
        {
            namespace fn = test_fn::func_2d;
            return std::make_tuple(InitFunction<2>{fn::density},
                                   InitFunctionArray{fn::vx, fn::vy, fn::vz},
                                   InitFunctionArray{fn::vthx, fn::vthy, fn::vthz});
        }
    }

    // loads every patch of the tiling and returns the particles sorted by cell, in load order
    static CellParticles load(std::vector<Box<int, dim>> const& tiling,
                              std::optional<std::size_t> seed, std::string const& pop = "protons",
                              int const level = 0, std::uint32_t const ppc = 20)
    {
        auto const [n, v, vth] = functions();
        MaxwellianParticleInitializer<ParticleArrayT, GridLayoutT> initializer{
            n, v, vth, ppc, seed, Basis::Cartesian, {nullptr, nullptr, nullptr}, 1e-5, pop};

        auto const dl = for_N_make_array<dim>([](auto) { return 0.1; });
        CellParticles cells;
        for (auto const& box : tiling)
        {
            auto const nbrCells = for_N_make_array<dim>([&](auto i) {
                return static_cast<std::uint32_t>(box.upper[i] - box.lower[i] + 1);
            });
            auto const origin = for_N_make_array<dim>([&](auto i) { return box.lower[i] * dl[i]; });
            GridLayoutT layout{dl, nbrCells, Point<double, dim>{origin}, box, level};
            ParticleArrayT particles{layout.AMRBox()};
            initializer.loadParticles(particles, layout);
            for (auto const& particle : particles)
                cells[particle.iCell].push_back(particle);
        }
        return cells;
    }
};

using Loader1D = SeededLoader<1>;
using Loader2D = SeededLoader<2>;
using Box1D    = Box<int, 1>;
using Box2D    = Box<int, 2>;

auto const wholeDomain1D = std::vector{Box1D{{0}, {99}}};


TEST(PerCellSeed, sameCellOfTwoPatchesGetsDifferentParticles)
{
    // regression: with a seed, every patch used to restart the same random stream
    auto const cells = Loader1D::load({Box1D{{0}, {19}}, Box1D{{20}, {39}}}, 1337);
    auto const& a    = cells.at({0});
    auto const& b    = cells.at({20});
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        EXPECT_NE(a[i].delta, b[i].delta);
}


TEST(PerCellSeed, particlesDoNotDependOnTheTiling1D)
{
    auto const reference = Loader1D::load(wholeDomain1D, 1337);
    ASSERT_EQ(reference.size(), 100u);

    for (auto const& tiling : {std::vector{Box1D{{0}, {49}}, Box1D{{50}, {99}}},
                               std::vector{Box1D{{37}, {99}}, Box1D{{0}, {36}}}})
        EXPECT_EQ(reference, Loader1D::load(tiling, 1337));
}


TEST(PerCellSeed, particlesDoNotDependOnTheTiling2D)
{
    auto const reference = Loader2D::load({Box2D{{0, 0}, {9, 9}}}, 1337);
    ASSERT_EQ(reference.size(), 100u);

    auto const tiled = Loader2D::load({Box2D{{0, 0}, {3, 6}}, Box2D{{4, 0}, {9, 6}},
                                       Box2D{{0, 7}, {6, 9}}, Box2D{{7, 7}, {9, 9}}},
                                      1337);
    EXPECT_EQ(reference, tiled);
}


TEST(PerCellSeed, neighbourCellsGetDifferentParticles)
{
    auto const cells = Loader2D::load({Box2D{{0, 0}, {9, 9}}}, 1337);
    for (auto const& [cell, particles] : cells)
        for (auto const& neighbour :
             {std::array{cell[0] + 1, cell[1]}, std::array{cell[0], cell[1] + 1}})
        {
            if (cells.count(neighbour))
            {
                EXPECT_NE(particles[0].delta, cells.at(neighbour)[0].delta);
            }
        }
}


TEST(PerCellSeed, sameSeedGivesSameParticlesOtherSeedDoesNot)
{
    auto const reference = Loader1D::load(wholeDomain1D, 1337);
    EXPECT_EQ(reference, Loader1D::load(wholeDomain1D, 1337));
    EXPECT_NE(reference.at({5})[0].delta, Loader1D::load(wholeDomain1D, 1338).at({5})[0].delta);
}


TEST(PerCellSeed, populationAndLevelChangeTheParticles)
{
    auto const reference = Loader1D::load(wholeDomain1D, 1337, "protons", 0);
    EXPECT_NE(reference.at({5})[0].delta,
              Loader1D::load(wholeDomain1D, 1337, "alpha", 0).at({5})[0].delta);
    EXPECT_NE(reference.at({5})[0].delta,
              Loader1D::load(wholeDomain1D, 1337, "protons", 1).at({5})[0].delta);
}


TEST(PerCellSeed, unseededLoadsDiffer)
{
    auto const a = Loader1D::load(wholeDomain1D, std::nullopt);
    auto const b = Loader1D::load(wholeDomain1D, std::nullopt);
    EXPECT_NE(a.at({5})[0].delta, b.at({5})[0].delta);
}


TEST(PerCellSeed, velocitiesAreMaxwellianAndUncorrelatedBetweenCells)
{
    // test functions give V = Vth = x at the cell center, so (v - V) / Vth ~ N(0, 1)
    std::uint32_t constexpr ppc = 1000;
    auto const cells            = Loader1D::load(wholeDomain1D, 1337, "protons", 0, ppc);

    double sum = 0, sum2 = 0, nbr = 0;
    std::vector<double> firstOfCell;
    for (auto const& [cell, particles] : cells)
    {
        double const x = (cell[0] + 0.5) * 0.1;
        for (auto const& particle : particles)
            for (auto const v : particle.v)
            {
                double const z = (v - x) / x;
                sum += z;
                sum2 += z * z;
                ++nbr;
            }
        firstOfCell.push_back((particles[0].v[0] - x) / x);
    }
    double const mean = sum / nbr;
    EXPECT_LT(std::abs(mean), 5. / std::sqrt(nbr));
    EXPECT_NEAR(std::sqrt(sum2 / nbr - mean * mean), 1., 0.05);

    double correlation = 0;
    for (std::size_t i = 0; i + 1 < firstOfCell.size(); ++i)
        correlation += firstOfCell[i] * firstOfCell[i + 1];
    correlation /= firstOfCell.size() - 1;
    EXPECT_LT(std::abs(correlation), 5. / std::sqrt(firstOfCell.size()));
}


TEST(PerCellSeed, keyDerivationIsPinned)
{
    // expected values computed independently (python), changing them changes every seeded run
    EXPECT_EQ(fnv1a64("protons"), 0xc59858089760c7faULL);
    EXPECT_EQ(cell_key(1337, 0xc59858089760c7faULL, 0, std::array{3}), 0x474d24db538d4c52ULL);
    EXPECT_EQ(cell_key(1337, 0xc59858089760c7faULL, 0, std::array{-2}), 0x47286a4d6bfd0c63ULL);
    EXPECT_EQ(cell_key(1337, 0xc59858089760c7faULL, 0, std::array{3, 7}), 0x16a13c97992bee23ULL);

    SplitMix64 generator{0x474d24db538d4c52ULL};
    for (auto const expected : {0x91d2d8d8dc5315bdULL, 0x9ca2dacc516034f6ULL, 0xca3deab3a2df2abaULL,
                                0xe78b14ebeeaa4a56ULL})
        EXPECT_EQ(generator(), expected);

    // reference vector of the SplitMix64 algorithm (seed 1234567)
    SplitMix64 reference{1234567};
    for (auto const expected :
         {6457827717110365317ULL, 3203168211198807973ULL, 9817491932198370423ULL})
        EXPECT_EQ(reference(), expected);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
