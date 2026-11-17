// The fully kinetic J and rho on a refined level (P4 plan D5): the level ghost band holds the
// deposit of the level ghost particles, made after the border sum, and the same-level ghost fills
// of the filter never overwrite it.

#include "phare_mpi.hpp"

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/fk_electric.hpp"
#include "core/data/field/field_traits.hpp"

#include "tests/simulator/per_test.hpp"

#include <set>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>


using namespace PHARE;

static constexpr std::size_t interp = PHARE_TEST_INTERP_ORDER;
static constexpr SimOpts opts{
    .dimension = 2, .interp_order = interp, .nbRefinedPart = 4, .fully_kinetic = true};
using Simulator_t = SimulatorTestParam<opts>;
using FKModel     = solver::PHARE_Types<opts>::FK::Model_t;
using GridLayout  = FKModel::gridlayout_type;
using Key         = std::array<int, 2>;

static constexpr int fineLevel = 1;
static constexpr double dtFine = 0.002 / 2; // job time step over the refinement ratio

// one simulator per process, destroyed before SAMRAI is shut down (main)
std::unique_ptr<Simulator_t> sim_;
Simulator_t& simulator()
{
    if (!sim_)
        sim_ = std::make_unique<Simulator_t>("job_o" + std::to_string(interp));
    return *sim_;
}

auto& level()
{
    auto& sim = simulator();
    if (sim.hierarchy->getNumberOfLevels() != 2)
        throw std::runtime_error("the job must give two levels");
    return *sim.hierarchy->getPatchLevel(fineLevel);
}

auto& messenger()
{
    return dynamic_cast<amr::HybridMessenger<FKModel>&>(
        simulator().getMultiPhysicsIntegrator()->messenger(fineLevel));
}


// the components of a vector field, or the field itself
template<typename Resource>
auto components(Resource& resource)
{
    if constexpr (core::IsField<Resource>)
        return std::array{&resource};
    else
        return std::array{&resource[0], &resource[1], &resource[2]};
}

Key amrIndex(GridLayout const& layout, auto i, auto j)
{
    auto const amr = layout.localToAMR(
        core::Point{static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j)}.as_signed());
    return {amr[0], amr[1]};
}

// local index of the node in [physical start - layers, physical end + layers] in both directions
bool within(GridLayout const& layout, auto const& field, std::uint32_t const layers, auto i, auto j)
{
    std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(i),
                                             static_cast<std::uint32_t>(j)};
    for (std::size_t d = 0; d < 2; ++d)
    {
        auto const direction = static_cast<core::Direction>(d);
        if (index[d] + layers < layout.physicalStartIndex(field, direction)
            or index[d] > layout.physicalEndIndex(field, direction) + layers)
            return false;
    }
    return true;
}



// T3a. Physical nodes hold a function, every ghost node NaN. After the ghost fill, the ghost
// nodes owned by another fine patch hold the function (corners included) and the level ghost band
// (no fine patch owns it) is still NaN: the fill is same-level only. The band must be NaN here:
// the field refine operators only write NaN nodes (composite_field_refiner.hpp), so a refine from
// the coarser level is only visible on NaN band nodes.
template<typename Resource, typename FillGhosts>
void checkGhostFillKeepsTheBand(Resource& resource, FillGhosts&& fillGhosts)
{
    auto& rm       = *simulator().getFKModel()->resourcesManager;
    auto& lvl      = level();
    auto const nan = std::numeric_limits<double>::quiet_NaN();

    auto const f = [](std::size_t c, auto const& r) {
        return std::sin(7.1 * r[0] + c) * std::cos(5.3 * r[1] - 0.3 * c) + c;
    };
    auto const position = [](auto const& layout, auto const& field, auto i, auto j) {
        return layout.fieldNodeCoordinates(
            field, layout.localToAMR(
                       core::Point{static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j)}
                           .as_signed()));
    };

    std::array<std::set<Key>, core::IsField<Resource> ? 1 : 3> owned;
    for (auto& patch : rm.enumerate(lvl, resource))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        std::size_t c     = 0;
        for (auto* field : components(resource))
        {
            layout.evalOnGhostBox(*field, [&](auto i, auto j) {
                if (within(layout, *field, 0, i, j))
                    owned[c].insert(amrIndex(layout, i, j));
            });
            ++c;
        }
    }

    for (auto& patch : rm.enumerate(lvl, resource))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        std::size_t c     = 0;
        for (auto* field : components(resource))
        {
            layout.evalOnGhostBox(*field, [&](auto i, auto j) {
                (*field)(i, j)
                    = within(layout, *field, 0, i, j) ? f(c, position(layout, *field, i, j)) : nan;
            });
            ++c;
        }
    }

    fillGhosts(lvl);

    std::size_t ownedGhosts = 0, band = 0;
    for (auto& patch : rm.enumerate(lvl, resource))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        std::size_t c     = 0;
        for (auto* field : components(resource))
        {
            layout.evalOnGhostBox(*field, [&](auto i, auto j) {
                auto const value = (*field)(i, j);
                if (owned[c].contains(amrIndex(layout, i, j)))
                {
                    EXPECT_NEAR(value, f(c, position(layout, *field, i, j)), 1e-14)
                        << field->name();
                    ownedGhosts += !within(layout, *field, 0, i, j);
                }
                else
                {
                    EXPECT_TRUE(std::isnan(value)) << field->name() << " level ghost band";
                    ++band;
                }
            });
            ++c;
        }
    }
    EXPECT_GT(ownedGhosts, 0u);
    EXPECT_GT(band, 0u);
}


TEST(FullyKineticLevelGhosts, filterGhostFillsAreSameLevelOnlyAndKeepTheLevelGhostBand)
{
    auto& state = simulator().getFKModel()->state;
    checkGhostFillKeepsTheBand(state.J,
                               [&](auto& lvl) { messenger().fillCurrentGhosts(state.J, lvl, 0.); });
    checkGhostFillKeepsTheBand(
        state.rho, [&](auto& lvl) { messenger().fillChargeDensityGhosts(state.rho, lvl, 0.); });
}



// the B-spline of the interpolation order, centered on the node, written here independently of
// the code under test
double shape(double x)
{
    x = std::abs(x);
    if constexpr (interp == 1)
        return x < 1 ? 1 - x : 0;
    else if constexpr (interp == 2)
        return x < 0.5 ? 0.75 - x * x : x < 1.5 ? 0.5 * (1.5 - x) * (1.5 - x) : 0;
    else
        return x < 1   ? 2. / 3. - x * x + 0.5 * x * x * x
               : x < 2 ? (2 - x) * (2 - x) * (2 - x) / 6.
                       : 0;
}


// T3b. J and rho deposited on the fine level by depositFKCurrent (no filter) against a direct sum
// over the particles each patch should see: all the fine domain particles (the border sum brings
// those of the neighbors), and the patch's own copies of the level ghost particles, each once.
//  - rho at x - dt/2 v on the physical nodes and the first ghost layer (the filter reads it)
//  - charge conservation on the physical nodes, rho(x + dt/2 v) - rho(x - dt/2 v) + dt div J = 0,
//    which needs J complete on them and on the first dual ghost layer
TEST(FullyKineticLevelGhosts, levelGhostParticlesDepositOnceAfterTheBorderSum)
{
    auto& model = *simulator().getFKModel();
    auto& rm    = *model.resourcesManager;
    auto& lvl   = level();
    auto& ions  = model.state.ions;
    auto& J     = model.state.J;
    auto& rho   = model.state.rho;
    ASSERT_GT(model.divergence_cleaning, 0.);
    ASSERT_EQ(model.current_filter_passes, 0u);

    solver::depositFKCurrent(lvl, model, messenger(), dtFine, 1., 1., 0.);

    struct Charge
    {
        std::array<int, 2> iCell;
        std::array<double, 2> delta;
        std::array<double, 2> v;
        double qw;
    };
    auto const charges = [](auto const& pop, auto const& particles, auto& out) {
        for (auto const& p : particles)
            out.push_back({p.iCell,
                           {static_cast<double>(p.delta[0]), static_cast<double>(p.delta[1])},
                           {p.v[0], p.v[1]},
                           pop.charge() * p.weight});
    };

    std::vector<Charge> domain;
    std::vector<std::vector<Charge>> levelGhosts;
    std::size_t nbrLevelGhosts = 0;
    for (auto& patch : rm.enumerate(lvl, ions))
    {
        levelGhosts.emplace_back();
        for (auto const& pop : ions)
        {
            charges(pop, pop.domainParticles(), domain);
            charges(pop, pop.levelGhostParticles(), levelGhosts.back());
        }
        nbrLevelGhosts += levelGhosts.back().size();
    }
    ASSERT_GT(nbrLevelGhosts, 0u);

    double constexpr halfSupport = 0.5 * (interp + 1);
    double maxRho = 0, maxRhoError = 0, maxResidual = 0;
    std::size_t nbrRhoNodes = 0, nbrConservationNodes = 0, patchIndex = 0;

    for (auto& patch : rm.enumerate(lvl, J, rho))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto const& dl    = layout.meshSize();
        auto const shape_ = rho.shape();
        std::vector<double> start(shape_[0] * shape_[1]), end(start.size());

        auto const deposit = [&](Charge const& p, double const timeShift, auto& out) {
            std::array<double, 2> pos; // in local node units, primal node i at i
            for (std::size_t d = 0; d < 2; ++d)
                pos[d] = p.iCell[d] - layout.AMRBox().lower[d]
                         + static_cast<int>(layout.physicalStartIndex(
                             core::QtyCentering::primal, static_cast<core::Direction>(d)))
                         + p.delta[d] + timeShift * p.v[d] / dl[d];
            auto const lo = [&](std::size_t d) {
                return std::max(0, static_cast<int>(std::ceil(pos[d] - halfSupport)));
            };
            auto const hi = [&](std::size_t d) {
                return std::min(static_cast<int>(shape_[d]) - 1,
                                static_cast<int>(std::floor(pos[d] + halfSupport)));
            };
            for (int i = lo(0); i <= hi(0); ++i)
                for (int j = lo(1); j <= hi(1); ++j)
                    out[i * shape_[1] + j] += p.qw * shape(pos[0] - i) * shape(pos[1] - j);
        };
        for (auto const* set : {&domain, &levelGhosts[patchIndex]})
            for (auto const& p : *set)
            {
                deposit(p, -0.5 * dtFine, start);
                deposit(p, +0.5 * dtFine, end);
            }

        layout.evalOnGhostBox(rho, [&](auto i, auto j) {
            if (!within(layout, rho, 1, i, j))
                return;
            auto const expected = start[i * shape_[1] + j];
            maxRho              = std::max(maxRho, std::abs(expected));
            maxRhoError         = std::max(maxRhoError, std::abs(rho(i, j) - expected));
            ++nbrRhoNodes;
        });

        layout.evalOnBox(rho, [&](auto i, auto j) {
            auto const k        = i * shape_[1] + j;
            auto const residual = end[k] - start[k] + dtFine * (J[0](i, j) - J[0](i - 1, j)) / dl[0]
                                  + dtFine * (J[1](i, j) - J[1](i, j - 1)) / dl[1];
            maxResidual = std::max(maxResidual, std::abs(residual));
            ++nbrConservationNodes;
        });
        ++patchIndex;
    }

    std::cout << "order " << interp << ": " << domain.size() << " domain and " << nbrLevelGhosts
              << " level ghost particles, rho on " << nbrRhoNodes << " nodes: max error "
              << maxRhoError / maxRho << " of max rho, charge conservation on "
              << nbrConservationNodes << " nodes: max residual " << maxResidual / maxRho << "\n";
    EXPECT_LT(maxRhoError, 1e-13 * maxRho);
    EXPECT_LT(maxResidual, 1e-13 * maxRho);
}


// T7a. E time labels (P4 plan D1, D8). The coarse E holds a + b tau, a and b linear in x and y, at
// three times: Eold_ (the messenger's copy, prepareStep) at t - Dt/2, the model E at t + Dt/2, the
// E one step ahead (the solver's) at t + 3 Dt/2. Eold_ must take E's own label. The fine level
// ghost band filled at t - Dt/4, t + Dt/4 (between Eold_ and E) and t + 3 Dt/4 (between E and the
// ahead E) must then hold a + b tau exactly: the linear refinement and time interpolation are exact
// on it.
TEST(FullyKineticLevelGhosts, fineElectricGhostsInterpolateBetweenTheThreeCoarseElectricFields)
{
    using VecField = FKModel::vecfield_type;
    auto& model    = *simulator().getFKModel();
    auto& rm       = *model.resourcesManager;
    auto& coarse   = *simulator().hierarchy->getPatchLevel(0);
    auto& fine     = level();
    auto& E        = model.state.electromag.E;

    // the messenger's and the solver's fields, by name
    VecField Eold{"FullyKineticModel-FullyKineticModel_Eold", core::HybridQuantity::Vector::E};
    VecField Eahead{"FK_Eahead", core::HybridQuantity::Vector::E};

    double constexpr t = 0., Dt = 0.002;
    auto const value = [](std::size_t c, auto const& r, double const tau) {
        auto const a = 1. + c + 2.1 * r[0] - 1.3 * r[1];
        auto const b = 400. * (1. + c + 0.7 * r[0] + 0.4 * r[1]);
        return a + b * tau;
    };
    auto const set = [&](auto& field, double const tau) {
        for (auto& patch : rm.enumerate(coarse, field))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            for (std::size_t c = 0; c < 3; ++c)
                layout.evalOnGhostBox(field[c], [&](auto i, auto j) {
                    field[c](i, j) = value(
                        c,
                        layout.fieldNodeCoordinates(
                            field[c], layout.localToAMR(core::Point{static_cast<std::uint32_t>(i),
                                                                    static_cast<std::uint32_t>(j)}
                                                            .as_signed())),
                        tau);
                });
            rm.setTime(field, *patch, tau);
        }
    };

    // the band: fine ghost nodes no fine patch owns
    std::array<std::set<Key>, 3> owned;
    for (auto& patch : rm.enumerate(fine, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        for (std::size_t c = 0; c < 3; ++c)
            layout.evalOnGhostBox(E[c], [&](auto i, auto j) {
                if (within(layout, E[c], 0, i, j))
                    owned[c].insert(amrIndex(layout, i, j));
            });
    }
    auto const checkBand = [&](double const tau) {
        messenger().fillElectricGhosts(E, fine, tau);
        double maxError = 0, scale = 0;
        std::size_t band = 0;
        for (auto& patch : rm.enumerate(fine, E))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            for (std::size_t c = 0; c < 3; ++c)
                layout.evalOnGhostBox(E[c], [&](auto i, auto j) {
                    if (owned[c].contains(amrIndex(layout, i, j)))
                        return;
                    auto const expected = value(
                        c,
                        layout.fieldNodeCoordinates(
                            E[c], layout.localToAMR(core::Point{static_cast<std::uint32_t>(i),
                                                                static_cast<std::uint32_t>(j)}
                                                        .as_signed())),
                        tau);
                    maxError = std::max(maxError, std::abs(E[c](i, j) - expected));
                    scale    = std::max(scale, std::abs(expected));
                    ++band;
                });
        }
        std::cout << "fine E band filled at t + " << (tau - t) / Dt << " Dt: " << band
                  << " nodes, max error " << maxError / scale << " of max |E|" << std::endl;
        EXPECT_GT(band, 0u);
        EXPECT_LT(maxError, 1e-12 * scale);
    };

    set(E, t - 0.5 * Dt);
    messenger().prepareStep(model, coarse, t);
    for (auto& patch : rm.enumerate(coarse, Eold))
        for (auto const time : rm.getTimes(Eold, *patch))
            EXPECT_EQ(time, t - 0.5 * Dt); // E's own label (D8), not currentTime

    set(E, t + 0.5 * Dt);
    messenger().firstStep(model, fine, t, t, t + Dt);
    checkBand(t - 0.25 * Dt);
    checkBand(t + 0.25 * Dt);

    set(Eahead, t + 1.5 * Dt);
    checkBand(t + 0.75 * Dt);
}


// T3e. The nodes the cleaning skips with boundary_cleaning_mask: the primal nodes of a patch whose
// 2^dim cells are some in a region and some not. The region is the level on a refined level (the
// patch and its same-level neighbours), the cells the finer level covers on a coarser one.
using Box2 = core::Box<int, 2>;

std::set<Key> nodeSet(std::vector<Box2> const& boxes)
{
    std::set<Key> nodes;
    for (auto const& box : boxes)
        core::for_each_point(box, [&](auto const& node) { nodes.insert({node[0], node[1]}); });
    return nodes;
}

// every node of the patch, one by one (node i is the lower corner of cell i)
template<typename InRegion>
std::set<Key> mixedNodes(Box2 const& cells, InRegion&& inRegion)
{
    std::set<Key> nodes;
    for (int i = cells.lower[0]; i <= cells.upper[0] + 1; ++i)
        for (int j = cells.lower[1]; j <= cells.upper[1] + 1; ++j)
        {
            int in = 0;
            for (int const di : {-1, 0})
                for (int const dj : {-1, 0})
                    in += inRegion(core::Point<int, 2>{i + di, j + dj});
            if (in > 0 and in < 4)
                nodes.insert({i, j});
        }
    return nodes;
}

Box2 box(int x0, int y0, int x1, int y1)
{
    return {core::Point<int, 2>{x0, y0}, core::Point<int, 2>{x1, y1}};
}

TEST(FullyKineticCleaningMask, regionEdgeNodesHaveCellsInAndOutOfTheRegion)
{
    auto const check = [](Box2 const& cells, std::vector<Box2> const& region) {
        auto const nodes = nodeSet(solver::regionEdgeNodes(cells, region));
        EXPECT_EQ(nodes, mixedNodes(cells, [&](auto const& c) { return core::isIn(c, region); }));
        return nodes;
    };

    // L-shaped level: the patch [0,7]^2 and a same-level neighbour [8,15]x[0,3], seen through the
    // particle ghost band. The face x = 8 is on the level boundary from y = 4 on only
    {
        auto const nodes = check(box(0, 0, 7, 7), {box(0, 0, 7, 7), box(8, 0, 9, 3)});
        EXPECT_EQ(nodes.size(), 32u - 3u);
        for (int const j : {1, 2, 3})
            EXPECT_FALSE(nodes.contains({8, j}));
        for (int const j : {0, 4, 8}) // (8, 4): the inner corner, on the boundary too
            EXPECT_TRUE(nodes.contains({8, j}));
    }
    // periodic in x over 16 cells: the neighbour [8,15]x[0,7] beyond the upper x face, its
    // periodic image [-8,-1]x[0,7] beyond the lower one. Only the faces normal to y are left
    {
        auto const nodes
            = check(box(0, 0, 7, 7), {box(0, 0, 7, 7), box(8, 0, 9, 7), box(-2, 0, -1, 7)});
        EXPECT_EQ(nodes.size(), 18u);
        for (auto const& [i, j] : nodes)
            EXPECT_TRUE(j == 0 or j == 8) << i << " " << j;
    }
    // a patch inside its level: nothing
    EXPECT_TRUE(check(box(0, 0, 7, 7), {box(0, 0, 7, 7), box(-2, -2, 9, 9)}).empty());

    // a coarser level patch [0,15]^2: two finer boxes sharing the face x = 8 for y in [4, 5], a
    // third one across the patch's upper x face (in the next patch too)
    {
        auto const nodes = check(box(0, 0, 15, 15),
                                 {box(4, 4, 7, 7), box(8, 4, 11, 5), box(14, 10, 17, 11)});
        EXPECT_FALSE(nodes.contains({8, 5})); // between two finer boxes
        EXPECT_TRUE(nodes.contains({8, 6}));
        EXPECT_FALSE(nodes.contains({6, 6})); // covered
        EXPECT_TRUE(nodes.contains({16, 10}));
        EXPECT_FALSE(nodes.contains({16, 11})); // the next patch is covered there too
        for (auto const& [i, j] : nodes)
            EXPECT_TRUE(i >= 0 and i <= 16 and j >= 0 and j <= 16) << i << " " << j;
    }
    EXPECT_TRUE(check(box(0, 0, 15, 15), {}).empty());
}


// the same on the job's levels, from the boxes SAMRAI gives: the root is 16^2 cells in 8^2
// patches, the fine level covers the root cells [4,11]^2, the fine cells [8,23]^2 in 4 patches.
// One process here: boxes of other ranks come through the same globalized box level, the
// periodic images one period away are added.
TEST(FullyKineticCleaningMask, skippedNodesOnTheJobLevelsAreTheEdgesOfTheFineLevel)
{
    auto& hierarchy = *simulator().hierarchy;
    level(); // two levels

    auto const check = [&](int const ln, auto const& skipped, Box2 const& region) {
        std::set<Key> all;
        for (auto const& patch : *hierarchy.getPatchLevel(ln))
        {
            auto const id    = amr::to_string(patch->getGlobalId());
            auto const cells = amr::phare_box_from<2>(patch->getBox());
            auto const nodes = skipped.contains(id) ? nodeSet(skipped.at(id)) : std::set<Key>{};
            EXPECT_EQ(nodes, mixedNodes(cells, [&](auto const& c) { return core::isIn(c, region); }))
                << "level " << ln << " patch " << cells;
            all.insert(nodes.begin(), nodes.end());
        }
        return all.size();
    };

    // fine level: its edge, from the same-level neighbours
    auto const fine = solver::cleaningSkippedNodes<GridLayout>(hierarchy, fineLevel, {});
    EXPECT_EQ(check(fineLevel, fine, box(8, 8, 23, 23)), 4u * 16u);

    // root: no level edge, the edge of the cells the fine level covers
    auto const covered = solver::finerLevelCoverage<2>(hierarchy, 0);
    auto const root    = solver::cleaningSkippedNodes<GridLayout>(hierarchy, 0, covered);
    EXPECT_EQ(check(0, root, box(4, 4, 11, 11)), 4u * 8u);

    // the finest level covers nothing
    EXPECT_TRUE(solver::finerLevelCoverage<2>(hierarchy, fineLevel).empty());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    auto const result = RUN_ALL_TESTS();
    sim_.reset();
    return result;
}
