// P5a: a fully kinetic level regridded during the run. One simulator per process: T14 runs first on
// the physical state of the run, T13 then overwrites the fields with test values.
//
// T14. The Gauss gap of a regridded level and the cleaning passes. After a few coarse steps the
// fine level is moved by two coarse cells through the integrator, as SAMRAI does on a regrid, with
// no cleaning pass at regrid; then F = div E - c^2 rho is measured on the nodes of the new level
// at least MARGIN fine cells inside its boundary, before and after each Marder pass (s = 1/2, rho
// fixed: cleanElectric, as the level initializer). On this square mesh a pass multiplies the mode
// (a, b) = (k_x dx, k_y dy) by 1 - (sin^2(a/2) + sin^2(b/2)) / 2. Checked:
// - every pass, node by node: E after the pass is E before it plus d dt grad F, F = div E - c^2 rho
//   set to 0 on the skipped nodes (core::Marder, done here by hand on the E the pass reads), on
//   the physical nodes of E_x, E_y, E_z unchanged, to rounding (1e-13 of max |E|);
// - a regression number: after K passes, the grid scale part F - S F (S one binomial pass
//   (1/4, 1/2, 1/4) along each direction) is at most GRID_LEFT of its value at regrid. It is not
//   G^K: F - S F also holds medium wavelengths, which a pass damps by about 0.9 (measured
//   0.338, 0.352, 0.366 at orders 1, 2, 3).
// The smooth part S F is reported: 10 passes leave about 62 percent of it (0.620, 0.638, 0.656),
// which the run-time passes, one per step, then keep cleaning. Marder passes at regrid do not
// close the smooth gap; a Poisson solve would (P4 plan pitfall 21).
//
// T13. The regridded E (P4 plan D6, P5a). The coarse E at t - Dt/2 and the coarse E one step
// ahead at t + Dt/2 are linear in x and y, the old fine level's E another linear function. A
// level regridded at t gets its E at t - Dt/4, its own half step (dt_f = Dt/2): on the physical
// nodes the old level had, the old level's E, copied (exact); on its other physical nodes and on
// its level ghosts, 3/4 E^{n-1/2} + 1/4 E^{n+1/2} refined, which the linear refinement gives to
// rounding (1e-14 of max |E|).

#include "phare_mpi.hpp"

#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/fk_electric.hpp"
#include "amr/level_initializer/fk_level_initializer.hpp"

#include "tests/simulator/per_test.hpp"

#include <SAMRAI/hier/Box.h>
#include <SAMRAI/hier/BoxLevel.h>
#include <SAMRAI/hier/BoxContainer.h>

#include <set>
#include <map>
#include <cmath>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>


using namespace PHARE;

static constexpr std::size_t interp = PHARE_TEST_INTERP_ORDER;
static constexpr SimOpts opts{
    .dimension = 2, .interp_order = interp, .nbRefinedPart = 4, .fully_kinetic = true};
using Simulator_t = SimulatorTestParam<opts>;
using FKModel     = solver::PHARE_Types<opts>::FK::Model_t;
using GridLayout  = FKModel::gridlayout_type;
using VecField    = FKModel::vecfield_type;
using Initializer = solver::FKLevelInitializer<FKModel>;
using Index2      = std::array<int, 2>;
using Key         = std::array<int, 3>; // component, AMR index

static constexpr int fineLevel    = 1;
static constexpr double Dt        = 0.002; // the job's time step, the fine one is Dt / 2
static constexpr std::size_t K    = 10;    // cleaning passes measured
static constexpr double GRID_LEFT = 0.40; // grid scale part left after K passes (measured <= 0.366)
static constexpr double STRENGTH  = 0.5;  // cleaning strength s of the passes
static constexpr int MARGIN       = 4;    // fine cells away from the level boundary
static constexpr std::size_t STEPS = 4;   // coarse steps before the regrid

// one simulator per process, destroyed before SAMRAI is shut down (main)
std::unique_ptr<Simulator_t> sim_;
Simulator_t& simulator()
{
    if (!sim_)
        sim_ = std::make_unique<Simulator_t>("job_regrid_o" + std::to_string(interp));
    return *sim_;
}

auto& messenger()
{
    return dynamic_cast<amr::HybridMessenger<FKModel>&>(
        simulator().getMultiPhysicsIntegrator()->messenger(fineLevel));
}

Index2 amrIndex(GridLayout const& layout, auto i, auto j)
{
    auto const amr = layout.localToAMR(
        core::Point{static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j)}.as_signed());
    return {amr[0], amr[1]};
}

auto position(GridLayout const& layout, auto const& field, auto i, auto j)
{
    return layout.fieldNodeCoordinates(
        field, layout.localToAMR(core::Point{static_cast<std::uint32_t>(i),
                                             static_cast<std::uint32_t>(j)}
                                     .as_signed()));
}

// a physical node of the field on this patch
bool physical(GridLayout const& layout, auto const& field, auto i, auto j)
{
    std::array<std::uint32_t, 2> const index{static_cast<std::uint32_t>(i),
                                             static_cast<std::uint32_t>(j)};
    for (std::size_t d = 0; d < 2; ++d)
    {
        auto const direction = static_cast<core::Direction>(d);
        if (index[d] < layout.physicalStartIndex(field, direction)
            or index[d] > layout.physicalEndIndex(field, direction))
            return false;
    }
    return true;
}

// replaces the fine level by one patch over the fine cells lower..upper, as SAMRAI does when the
// boxes of a level change (GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel), and
// returns the old level
std::shared_ptr<SAMRAI::hier::PatchLevel> replaceFineLevel(Index2 const lower, Index2 const upper)
{
    auto& hierarchy = *simulator().hierarchy;
    auto oldLevel   = hierarchy.getPatchLevel(fineLevel);

    SAMRAI::hier::Box const box{SAMRAI::hier::Index{std::vector<int>{lower[0], lower[1]}},
                                SAMRAI::hier::Index{std::vector<int>{upper[0], upper[1]}},
                                SAMRAI::hier::BlockId{0}, SAMRAI::hier::LocalId{0},
                                SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank()};
    SAMRAI::hier::BoxContainer boxes;
    boxes.pushBack(box);
    SAMRAI::hier::BoxLevel const boxLevel{boxes, oldLevel->getRatioToLevelZero(),
                                          oldLevel->getGridGeometry()};

    hierarchy.removePatchLevel(fineLevel);
    hierarchy.makeNewPatchLevel(fineLevel, boxLevel);

    // the new<==>coarser Connectors SAMRAI's regrid caches for the new level
    // (GriddingAlgorithm::regridFinerLevel_createAndInstallNewLevel): the regrid schedules, from
    // the old level not in the hierarchy, connect the new level to the coarser one through them
    auto const newLevel = hierarchy.getPatchLevel(fineLevel);
    newLevel->findConnectorWithTranspose(
        *hierarchy.getPatchLevel(fineLevel - 1),
        hierarchy.getRequiredConnectorWidth(fineLevel, fineLevel - 1, true),
        hierarchy.getRequiredConnectorWidth(fineLevel - 1, fineLevel, true),
        SAMRAI::hier::CONNECTOR_IMPLICIT_CREATION_RULE, true);
    return oldLevel;
}


// F = div E - c^2 rho on the primal nodes of the (one patch) fine level at least MARGIN cells
// inside its boundary
std::map<Index2, double> gaussResidual(FKModel& model)
{
    auto& rm    = *model.resourcesManager;
    auto& E     = model.state.electromag.E;
    auto& rho   = model.state.rho;
    auto& level = *simulator().hierarchy->getPatchLevel(fineLevel);

    std::map<Index2, double> F;
    for (auto& patch : rm.enumerate(level, E, rho))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto const& dl    = layout.meshSize();
        auto const& box   = patch->getBox();
        layout.evalOnBox(rho, [&](auto i, auto j) {
            auto const node = amrIndex(layout, i, j);
            for (std::size_t d = 0; d < 2; ++d)
                if (node[d] < box.lower()[d] + MARGIN or node[d] > box.upper()[d] + 1 - MARGIN)
                    return;
            // the dual node i of E_d is at i + 1/2 (core::Marder)
            F[node] = (E[0](i, j) - E[0](i - 1, j)) / dl[0]
                      + (E[1](i, j) - E[1](i, j - 1)) / dl[1] - model.c2() * rho(i, j);
        });
    }
    return F;
}

struct Parts
{
    double grid = 0, smooth = 0, total = 0; // 2-norms
    bool finite = true;
};

// grid scale part F - S F and smooth part S F, on the nodes whose 8 neighbors have F
Parts split(std::map<Index2, double> const& F)
{
    Parts parts;
    double constexpr w[3] = {0.25, 0.5, 0.25};
    for (auto const& [node, value] : F)
    {
        parts.finite = parts.finite and std::isfinite(value);
        parts.total += value * value;
        double smooth = 0;
        bool complete = true;
        for (int a = -1; a <= 1 and complete; ++a)
            for (int b = -1; b <= 1 and complete; ++b)
            {
                auto const it = F.find({node[0] + a, node[1] + b});
                complete      = it != F.end();
                if (complete)
                    smooth += w[a + 1] * w[b + 1] * it->second;
            }
        if (!complete)
            continue;
        parts.smooth += smooth * smooth;
        parts.grid += (value - smooth) * (value - smooth);
    }
    parts.grid   = std::sqrt(parts.grid);
    parts.smooth = std::sqrt(parts.smooth);
    parts.total  = std::sqrt(parts.total);
    return parts;
}


bool isSkipped(std::vector<core::Box<int, 2>> const& boxes, Index2 const& node)
{
    return std::any_of(boxes.begin(), boxes.end(), [&](auto const& box) {
        return box.lower[0] <= node[0] and node[0] <= box.upper[0] and box.lower[1] <= node[1]
               and node[1] <= box.upper[1];
    });
}

// one cleaning pass (cleanElectric, s = STRENGTH) on the fine level, checked node by node against
// one core::Marder step done here by hand on the E the pass reads (ghosts filled at `time`).
// Returns the largest |E - expected| over the physical nodes, relative to max |E|
double checkedPass(FKModel& model, VecField& E, double const dt, double const time,
                   solver::PatchBoxes<2> const& skipped)
{
    auto& rm    = *model.resourcesManager;
    auto& rho   = model.state.rho;
    auto& level = *simulator().hierarchy->getPatchLevel(fineLevel);

    messenger().fillElectricPatchGhosts(E, level, time);
    std::vector<std::array<std::vector<double>, 3>> expected; // per patch, per component
    for (auto& patch : rm.enumerate(level, E, rho))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto const& dl    = layout.meshSize();
        auto const ddt    = STRENGTH / (2. / (dl[0] * dl[0]) + 2. / (dl[1] * dl[1])); // d dt
        auto const it     = skipped.find(amr::to_string(patch->getGlobalId()));
        auto const ny     = rho.shape()[1];

        std::vector<double> F(rho.shape()[0] * ny, 0.);
        layout.evalOnBox(rho, [&](auto i, auto j) {
            if (it != skipped.end() and isSkipped(it->second, amrIndex(layout, i, j)))
                return;
            F[i * ny + j] = -model.c2() * rho(i, j) + (E[0](i, j) - E[0](i - 1, j)) / dl[0]
                            + (E[1](i, j) - E[1](i, j - 1)) / dl[1];
        });

        auto& out = expected.emplace_back();
        layout.evalOnBox(E[0], [&](auto i, auto j) {
            out[0].push_back(E[0](i, j) + ddt / dl[0] * (F[(i + 1) * ny + j] - F[i * ny + j]));
        });
        layout.evalOnBox(E[1], [&](auto i, auto j) {
            out[1].push_back(E[1](i, j) + ddt / dl[1] * (F[i * ny + j + 1] - F[i * ny + j]));
        });
        layout.evalOnBox(E[2], [&](auto i, auto j) { out[2].push_back(E[2](i, j)); });
    }

    solver::cleanElectric(level, model, messenger(), E, dt, STRENGTH, 1, time, skipped);

    double error = 0, scale = 0;
    std::size_t p = 0;
    for (auto& patch : rm.enumerate(level, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        for (std::size_t c = 0; c < 3; ++c)
        {
            std::size_t n = 0;
            layout.evalOnBox(E[c], [&](auto i, auto j) {
                error = std::max(error, std::abs(E[c](i, j) - expected[p][c][n]));
                scale = std::max(scale, std::abs(E[c](i, j)));
                ++n;
            });
            EXPECT_EQ(n, expected[p][c].size());
        }
        ++p;
    }
    return error / scale;
}


TEST(FullyKineticRegrid, gaussGapOfARegriddedLevelAndTheCleaningPasses)
{
    auto& sim        = simulator();
    auto& model      = *sim.getFKModel();
    auto& rm         = *model.resourcesManager;
    auto& integrator = *sim.getMultiPhysicsIntegrator();
    auto& E          = model.state.electromag.E;

    for (std::size_t n = 0; n < STEPS; ++n)
        sim.advance();
    auto const t       = sim.currentTime();
    auto const dt      = 0.5 * Dt;
    auto const regrids = integrator.regridCount(fineLevel);

    // the job's box covers fine cells 20..43 along both directions: moved by 2 coarse cells in x
    Initializer::regridCleaningPasses = 0;
    auto const oldLevel               = replaceFineLevel({24, 20}, {47, 43});
    integrator.initializeLevelData(sim.hierarchy, fineLevel, t, true, false, oldLevel);
    Initializer::regridCleaningPasses = 10;

    EXPECT_EQ(integrator.regridCount(fineLevel), regrids + 1);
    auto& level = *sim.hierarchy->getPatchLevel(fineLevel);
    for (auto& patch : rm.enumerate(level, E))
        for (auto const time : rm.getTimes(E, *patch))
            EXPECT_NEAR(time, t - 0.5 * dt, 1e-12) << "E at the level's own half step";

    // the nodes the initializer's passes skip (boundary_cleaning_mask: the level edge)
    auto const skipped
        = model.boundary_cleaning_mask
              ? solver::cleaningSkippedNodes<GridLayout>(*sim.hierarchy, fineLevel, {})
              : solver::PatchBoxes<GridLayout::dimension>{};
    std::vector<Parts> parts{split(gaussResidual(model))};
    std::vector<double> passErrors{0.};
    for (std::size_t k = 1; k <= K; ++k)
    {
        passErrors.push_back(checkedPass(model, E, dt, t - 0.5 * dt, skipped));
        parts.push_back(split(gaussResidual(model)));
    }

    std::cout << std::setprecision(3) << "order " << interp << ", regrid at t = " << t
              << ": |F| (total, grid scale, smooth) after k passes, relative to k = 0;"
              << " pass error relative to max |E|\n";
    for (std::size_t k = 0; k <= K; ++k)
    {
        ASSERT_TRUE(parts[k].finite) << "pass " << k;
        std::cout << "  k = " << k << ": " << parts[k].total / parts[0].total << " "
                  << parts[k].grid / parts[0].grid << " " << parts[k].smooth / parts[0].smooth
                  << "  (pass error " << passErrors[k] << ")" << std::endl;
        EXPECT_LE(passErrors[k], 1e-13) << "pass " << k;
    }
    std::cout << "  at k = 0: total " << parts[0].total << ", grid scale " << parts[0].grid
              << ", smooth " << parts[0].smooth << std::endl;
    ASSERT_GT(parts[0].grid, 0.);
    EXPECT_LE(parts[K].grid, GRID_LEFT * parts[0].grid);
}


TEST(FullyKineticRegrid, regriddedElectricIsCopiedWhereTheOldLevelWasInterpolatedElsewhere)
{
    auto& sim    = simulator();
    auto& model  = *sim.getFKModel();
    auto& rm     = *model.resourcesManager;
    auto& coarse = *sim.hierarchy->getPatchLevel(0);
    auto& E      = model.state.electromag.E;
    VecField Eahead{"FK_Eahead", core::HybridQuantity::Vector::E}; // the solver's, by name

    double constexpr t        = 0.;
    double constexpr fillTime = t - 0.25 * Dt;
    auto const coarseValue    = [](std::size_t c, auto const& r, double const tau) {
        auto const a = 1. + c + 2.1 * r[0] - 1.3 * r[1];
        auto const b = 400. * (1. + c + 0.7 * r[0] + 0.4 * r[1]);
        return a + b * tau;
    };
    auto const oldValue = [](std::size_t c, auto const& r) {
        return -7. - c + 0.9 * r[0] + 0.2 * r[1];
    };
    auto const set = [&](auto& level, auto& field, auto&& value, double const tau) {
        for (auto& patch : rm.enumerate(level, field))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            for (std::size_t c = 0; c < 3; ++c)
                layout.evalOnGhostBox(field[c], [&](auto i, auto j) {
                    field[c](i, j) = value(c, position(layout, field[c], i, j));
                });
            rm.setTime(field, *patch, tau);
        }
    };

    set(coarse, E, [&](auto c, auto const& r) { return coarseValue(c, r, t - 0.5 * Dt); },
        t - 0.5 * Dt);
    set(coarse, Eahead, [&](auto c, auto const& r) { return coarseValue(c, r, t + 0.5 * Dt); },
        t + 0.5 * Dt);

    // the old level (T14's: fine cells 24..47 x 20..43) and its physical nodes
    std::set<Key> oldNodes;
    auto& old = *sim.hierarchy->getPatchLevel(fineLevel);
    set(old, E, oldValue, fillTime);
    for (auto& patch : rm.enumerate(old, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        for (std::size_t c = 0; c < 3; ++c)
            layout.evalOnGhostBox(E[c], [&](auto i, auto j) {
                if (physical(layout, E[c], i, j))
                {
                    auto const node = amrIndex(layout, i, j);
                    oldNodes.insert({static_cast<int>(c), node[0], node[1]});
                }
            });
    }

    // moved by 4 coarse cells back in x: fine cells 16..23 are new, 24..39 were covered
    auto const oldLevel = replaceFineLevel({16, 20}, {39, 43});
    auto& fine          = *sim.hierarchy->getPatchLevel(fineLevel);
    for (auto& patch : fine)
        model.allocate(*patch, t);
    messenger().regridElectric(E, sim.hierarchy, fineLevel, oldLevel, fillTime);

    double copyError = 0, interpolationError = 0, scale = 0;
    std::size_t copied = 0, interpolated = 0, ghosts = 0;
    for (auto& patch : rm.enumerate(fine, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        for (std::size_t c = 0; c < 3; ++c)
            layout.evalOnGhostBox(E[c], [&](auto i, auto j) {
                auto const r    = position(layout, E[c], i, j);
                auto const node = amrIndex(layout, i, j);
                bool const isPhysical = physical(layout, E[c], i, j);
                if (isPhysical and oldNodes.contains({static_cast<int>(c), node[0], node[1]}))
                {
                    copyError = std::max(copyError, std::abs(E[c](i, j) - oldValue(c, r)));
                    ++copied;
                    return;
                }
                auto const expected = coarseValue(c, r, fillTime);
                interpolationError
                    = std::max(interpolationError, std::abs(E[c](i, j) - expected));
                scale = std::max(scale, std::abs(expected));
                ++(isPhysical ? interpolated : ghosts);
            });
    }

    std::cout << "order " << interp << ": regridded E, " << copied << " nodes copied (max error "
              << copyError << "), " << interpolated << " physical and " << ghosts
              << " ghost nodes interpolated (max error " << interpolationError / scale
              << " of max |E|)" << std::endl;
    EXPECT_GT(copied, 0u);
    EXPECT_GT(interpolated, 0u);
    EXPECT_GT(ghosts, 0u);
    EXPECT_EQ(copyError, 0.);
    EXPECT_LT(interpolationError, 1e-14 * scale);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    auto const result = RUN_ALL_TESTS();
    sim_.reset();
    return result;
}
