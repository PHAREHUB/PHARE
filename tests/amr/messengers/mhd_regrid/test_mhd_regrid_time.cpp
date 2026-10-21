#include "phare_mpi.hpp"

#include "tests/simulator/per_test.hpp"

#include <map>
#include <sstream>
#include <string>
#include <vector>

/*
 * A level regridded or created between two syncs of its coarser level must take its ghosts from
 * the coarser levels interpolated to the regrid time, as a ghost fill during a step does.
 *
 * After one level-0 step, level 1 holds its last substep's start and end, level 0 its step's.
 * Level 2 is regridded, keeping its boxes, or removed and created again, at the step end and at
 * two times within the last level-1 substep. Its ghosts must vary linearly with the regrid time:
 * they differ between the regrid times, and the regrid at the midpoint gives their average.
 */

using MHDSimulator = SimulatorTestParam<SimOpts{2, 0, 0, MHDOpts::ReconstructionType::WENOZ,
                                                MHDOpts::SlopeLimiterType::None,
                                                MHDOpts::RiemannSolverType::Rusanov}>;

// every value of B, rho, rhoV and Etot on a patch, ghosts included, keyed by the patch box
using LevelValues = std::map<std::string, std::vector<double>>;


LevelValues regridLevel2(MHDSimulator& sim, double const regridTime, bool const createLevel)
{
    auto& hierarchy = *sim.hierarchy;
    auto& model     = *sim.getMHDModel();
    auto& state     = model.state;

    if (createLevel)
        hierarchy.removePatchLevel(2);

    std::vector<int> const tagBuffer(hierarchy.getMaxNumberOfLevels(), 1);
    sim.getIntegrator()->griddingAlgorithm()->regridAllFinerLevels(1, tagBuffer, 0, regridTime);

    LevelValues values;
    for (auto& patch : *hierarchy.getPatchLevel(2))
    {
        auto onPatch = model.resourcesManager->setOnPatch(*patch, state.B, state.rho, state.rhoV,
                                                          state.Etot);
        std::ostringstream box;
        box << patch->getBox();
        auto& patchValues = values[box.str()];

        auto const append = [&](auto const& field) {
            patchValues.insert(patchValues.end(), field.data(), field.data() + field.size());
        };
        for (auto const& component : state.B)
            append(component);
        append(state.rho);
        for (auto const& component : state.rhoV)
            append(component);
        append(state.Etot);
    }
    return values;
}


double maxDifference(LevelValues const& values, LevelValues const& reference)
{
    EXPECT_EQ(values.size(), reference.size());
    double maxDiff = 0;
    for (auto const& [box, patchValues] : values)
    {
        auto const& referenceValues = reference.at(box);
        EXPECT_EQ(patchValues.size(), referenceValues.size());
        for (std::size_t i = 0; i < patchValues.size(); ++i)
            maxDiff = std::max(maxDiff, std::abs(patchValues[i] - referenceValues[i]));
    }
    return maxDiff;
}


LevelValues average(LevelValues a, LevelValues const& b)
{
    for (auto& [box, patchValues] : a)
        for (std::size_t i = 0; i < patchValues.size(); ++i)
            patchValues[i] = 0.5 * (patchValues[i] + b.at(box)[i]);
    return a;
}


void expectGhostsAtTheRegridTime(bool const createLevel)
{
    MHDSimulator sim;
    sim.advance();

    double const end = sim.currentTime();
    double const dt  = sim.timeStep();

    auto const atEnd   = regridLevel2(sim, end, createLevel);
    auto const before  = regridLevel2(sim, end - dt / 4, createLevel);
    auto const between = regridLevel2(sim, end - dt / 8, createLevel);

    ASSERT_FALSE(atEnd.empty());
    EXPECT_GT(maxDifference(before, atEnd), 1e-8);
    EXPECT_LT(maxDifference(between, average(before, atEnd)), 1e-12);
}


TEST(MHDRegrid, regriddedLevelTakesGhostsAtTheRegridTime)
{
    expectGhostsAtTheRegridTime(false);
}


TEST(MHDRegrid, createdLevelTakesGhostsAtTheRegridTime)
{
    expectGhostsAtTheRegridTime(true);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    return RUN_ALL_TESTS();
}
