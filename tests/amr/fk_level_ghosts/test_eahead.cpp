// T6 (P4 plan). The coarse level computes E one step ahead, Eahead = E^{n+3/2}, before its finer
// level's substeps: the newer end of the time interpolation of the fine E ghosts past E^{n+1/2}
// (D1). It is computed before the synchronization, which then changes E^{n+1/2} on the covered side
// (the fine time average is coarsened onto it) and B^{n+1} (reflux). This measures, on the coarse
// nodes the fine ghost refinement reads, how far Eahead is from the E^{n+3/2} the next coarse step
// computes after that synchronization, by row: one coarse cell inside the box, the interface, and
// every row outside the fine ghosts read. With the fine field ghosts 4, 6, 6 at orders 1 to 3
// (fe21e3bab), the outermost fine dual ghost is (g - 1/2)/2 = 1.75, 2.75, 2.75 coarse cells out and
// the composite refiner builds it from the coarse duals of its cell and of the next one out:
// 1.5 and 2.5 at order 1, 2.5 and 3.5 at orders 2 and 3 (the outermost primal ghost, g/2 = 2, 3, 3,
// is a copy of the coincident coarse node). So the rows read go up to 3 coarse cells outside at
// order 1 and 4 at orders 2 and 3: ceil(g/2 + 1/2), taken here from the fine layout.
// Only the rows outside measure the replay gap: no synchronization ever changes them. Inside and on
// the interface, the target E^{n+3/2} is taken after the next step, so it already includes that
// step's synchronization (the coarsened fine average), which the replay cannot know: there the
// difference mixes the two. Hence the sanity bound on the first row outside.
// It also compares, at the fill time t_n + 3/4 Dt, the fine ghost input of the scheme,
// 3/4 E^{n+1/2} + 1/4 Eahead, and of a zeroth-order hold, E^{n+1/2}, with the ideal one,
// 3/4 E^{n+1/2} + 1/4 E^{n+3/2}: errors 1/4 |Eahead - E^{n+3/2}| and 1/4 |E^{n+1/2} - E^{n+3/2}|,
// computed from the stored coarse fields, no switch in the code.
// The case is the job module on the command line (job_eahead_<case>.py): wave32 and wave16, a
// standing light wave in near vacuum (smooth part, second order: the ratio of the two should be
// about 4), thermal (noise part). Report only, but for one sanity bound on the smooth cases: the
// scheme beats the hold one row outside the box.

#include "phare_mpi.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include "tests/simulator/per_test.hpp"

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

static constexpr SimOpts opts{
    .dimension = 2, .interp_order = 1, .nbRefinedPart = 4, .fully_kinetic = true};
using Simulator_t = SimulatorTestParam<opts>;
using FKModel     = solver::PHARE_Types<opts>::FK::Model_t;
using GridLayout  = FKModel::gridlayout_type;
using VecField    = FKModel::vecfield_type;

static std::string job                = "job_eahead_wave32";
static constexpr std::size_t nbrSteps = 12;

// component and twice the position in coarse cells: the dual nodes at half integers
using Key = std::array<int, 3>;

// rows: 0 one cell inside, 1 the interface, 1 + k the k-th cell outside, k = 1 .. outsideRows
static constexpr int coveredRow = 0, interface = 1, outsideRow1 = 2;

std::string rowName(int const row)
{
    return row == coveredRow  ? "1 cell inside"
           : row == interface ? "interface"
                              : std::to_string(row - interface) + " cells outside";
}

std::unique_ptr<Simulator_t> sim_;


struct Box2
{
    std::array<double, 2> lo, hi; // the box faces, in coarse cell units
};

// the refined box on the coarse level, from the fine patches (ratio 2)
Box2 refinedBox(Simulator_t& sim)
{
    Box2 box{{1e300, 1e300}, {-1e300, -1e300}};
    for (auto const& patch : *sim.hierarchy->getPatchLevel(1))
        for (std::size_t d = 0; d < 2; ++d)
        {
            box.lo[d] = std::min(box.lo[d], patch->getBox().lower()[d] / 2.);
            box.hi[d] = std::max(box.hi[d], (patch->getBox().upper()[d] + 1) / 2.);
        }
    return box;
}

// the coarse rows outside the box the fine ghosts read, from the fine field ghost width g
int outsideRows(Simulator_t& sim)
{
    auto const& patch = *sim.hierarchy->getPatchLevel(1)->begin();
    auto const g      = amr::layoutFromPatch<GridLayout>(*patch).physicalStartIndex(
        core::QtyCentering::primal, core::Direction::X);
    return static_cast<int>(std::ceil(g / 2. + 0.5));
}

// the row of a node at position p (coarse cell units), -1 if not read by the fine ghosts
int rowOf(Box2 const& box, int const nbrOutside, std::array<double, 2> const& p)
{
    double outside = 0, inside = 1e300;
    for (std::size_t d = 0; d < 2; ++d)
    {
        outside = std::max({outside, box.lo[d] - p[d], p[d] - box.hi[d]});
        inside  = std::min({inside, p[d] - box.lo[d], box.hi[d] - p[d]});
    }
    if (outside > 0)
        return outside <= nbrOutside ? interface + static_cast<int>(std::ceil(outside)) : -1;
    if (inside == 0)
        return interface;
    return inside <= 1 ? coveredRow : -1;
}


struct Snapshot
{
    std::map<Key, double> values; // nodes in a row class
    std::map<Key, int> rows;
    double maxAbs = 0; // over all the coarse physical nodes
};

Snapshot snapshot(Simulator_t& sim, VecField& field, Box2 const& box, int const nbrOutside)
{
    auto& rm     = *sim.getFKModel()->resourcesManager;
    auto& coarse = *sim.hierarchy->getPatchLevel(0);
    Snapshot snap;
    for (auto& patch : rm.enumerate(coarse, field))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        auto const& dl    = layout.meshSize();
        for (std::size_t c = 0; c < 3; ++c)
            layout.evalOnBox(field[c], [&](auto i, auto j) {
                auto const value = field[c](i, j);
                snap.maxAbs      = std::max(snap.maxAbs, std::abs(value));
                auto const x     = layout.fieldNodeCoordinates(
                    field[c], layout.localToAMR(core::Point{static_cast<std::uint32_t>(i),
                                                            static_cast<std::uint32_t>(j)}
                                                        .as_signed()));
                std::array<double, 2> const p{x[0] / dl[0], x[1] / dl[1]};
                if (auto const row = rowOf(box, nbrOutside, p); row >= 0)
                {
                    Key const key{static_cast<int>(c), static_cast<int>(std::lround(2 * p[0])),
                                  static_cast<int>(std::lround(2 * p[1]))};
                    snap.values.emplace(key, value); // shared nodes: the first patch's
                    snap.rows.emplace(key, row);
                }
            });
    }
    return snap;
}


TEST(FullyKineticEahead, aheadElectricAgainstTheNextStepsSynchronizedElectric)
{
    sim_                  = std::make_unique<Simulator_t>(job);
    auto& sim             = *sim_;
    auto& E               = sim.getFKModel()->state.electromag.E;
    auto const box        = refinedBox(sim);
    auto const nbrOutside = outsideRows(sim);
    auto const nbrRows    = static_cast<std::size_t>(outsideRow1 + nbrOutside);
    VecField Eahead{"FK_Eahead", core::HybridQuantity::Vector::E}; // the solver's, by name

    // after step n: E^{n+1/2} synchronized, Eahead^{n+3/2} from before the synchronization
    std::vector<Snapshot> electric, ahead;
    for (std::size_t n = 0; n < nbrSteps; ++n)
    {
        sim.advance();
        electric.push_back(snapshot(sim, E, box, nbrOutside));
        ahead.push_back(snapshot(sim, Eahead, box, nbrOutside));
    }

    double scale = 0;
    for (auto const& snap : electric)
        scale = std::max(scale, snap.maxAbs);
    ASSERT_GT(scale, 0.);

    std::vector<double> schemeMax(nbrRows), schemeSq(nbrRows), holdSq(nbrRows);
    std::vector<std::size_t> count(nbrRows);
    for (std::size_t n = 0; n + 1 < nbrSteps; ++n)
        for (auto const& [key, row] : electric[n].rows)
        {
            auto const next   = electric[n + 1].values.at(key);
            auto const scheme = 0.25 * std::abs(ahead[n].values.at(key) - next);
            auto const hold   = 0.25 * std::abs(electric[n].values.at(key) - next);
            ASSERT_TRUE(std::isfinite(scheme) and std::isfinite(hold));
            schemeMax[row] = std::max(schemeMax[row], scheme);
            schemeSq[row] += scheme * scheme;
            holdSq[row] += hold * hold;
            ++count[row];
        }

    std::cout << std::setprecision(3) << job << ": fine ghost input error at t + 3/4 Dt, "
              << nbrSteps - 1 << " steps, relative to max |E| = " << scale << "\n";
    for (int row = 0; row < static_cast<int>(nbrRows); ++row)
    {
        ASSERT_GT(count[row], 0u) << rowName(row);
        auto const schemeRms = std::sqrt(schemeSq[row] / count[row]) / scale;
        auto const holdRms   = std::sqrt(holdSq[row] / count[row]) / scale;
        std::cout << "  " << rowName(row) << " (" << count[row] << " node-steps): scheme max "
                  << schemeMax[row] / scale << " rms " << schemeRms << ", hold rms " << holdRms
                  << ", hold / scheme " << holdRms / schemeRms << std::endl;
        if (job.find("wave") != std::string::npos and row == outsideRow1)
        {
            EXPECT_LT(schemeRms, holdRms); // sanity: the trick must matter on a smooth wave
        }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
        job = argv[1];
    PHARE::SamraiLifeCycle samsam(argc, argv);
    auto const result = RUN_ALL_TESTS();
    sim_.reset();
    return result;
}
