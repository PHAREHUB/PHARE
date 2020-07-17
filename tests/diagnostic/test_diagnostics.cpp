
// input/input_1d_ratio_2.txt is unused but a reference

#include <mpi.h>

#include "test_diagnostics.h"

using namespace PHARE;
using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;


template<typename Simulator, typename Hi5Diagnostic>
void validateFluidDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using namespace std::string_literals;

    using GridLayout  = typename Simulator::PHARETypes::GridLayout_t;
    auto& hybridModel = *sim.getHybridModel();

    auto checkF = [&](auto& layout, auto& path, auto tree, auto name, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name));
        checkField(hifile->file(), layout, val, path + name, FieldDomainFilter{});
    };
    auto checkVF = [&](auto& layout, auto& path, auto tree, auto name, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name));
        checkVecField(hifile->file(), layout, val, path + name, FieldDomainFilter{});
    };

    auto visit = [&](GridLayout& layout, std::string patchID, size_t iLevel) {
        auto path  = hi5.getPatchPath(iLevel, patchID);
        auto& ions = hi5.modelView.getIons();
        for (auto& pop : ions)
        {
            checkF(layout, path, "/ions/pop/" + pop.name(), "/density"s, pop.density());
            checkVF(layout, path, "/ions/pop/" + pop.name(), "/flux"s, pop.flux());
        }
        checkF(layout, path, "/ions"s, "/density"s, ions.density());

        std::string tree{"/ions"}, var{"/bulkVelocity"};
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + var));
        checkVecField(hifile->file(), layout, ions.velocity(), path + var, FieldDomainFilter{});
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


TYPED_TEST(SimulatorTest, fluid)
{
    TypeParam sim;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;

    { // Scoped to destruct after dump
        Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.fluid("/ions/density"))
            .addDiagDict(hi5.fluid("/ions/bulkVelocity"))
            .addDiagDict(hi5.fluid("/ions/pop/alpha/density"))
            .addDiagDict(hi5.fluid("/ions/pop/alpha/flux"))
            .addDiagDict(hi5.fluid("/ions/pop/protons/density"))
            .addDiagDict(hi5.fluid("/ions/pop/protons/flux"));
        sim.dump(hi5.dMan);
    }

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, HighFive::File::ReadOnly};
    validateFluidDump(sim, hi5);
}


template<typename Simulator, typename Hi5Diagnostic>
void validateElectromagDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();

    auto checkVF = [&](auto& layout, auto& path, auto tree, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree));
        checkVecField(hifile->file(), layout, val, path + tree);
    };

    auto visit = [&](GridLayout& layout, std::string patchID, size_t iLevel) {
        auto path = hi5.getPatchPath(iLevel, patchID) + "/";
        checkVF(layout, path, "/EM_B", hybridModel.state.electromag.B);
        checkVF(layout, path, "/EM_E", hybridModel.state.electromag.E);
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}

TYPED_TEST(SimulatorTest, electromag)
{
    TypeParam sim;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;
    { // scoped to destruct after dump
        Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.electromag("/EM_B")).addDiagDict(hi5.electromag("/EM_E"));
        sim.dump(hi5.dMan);
    }

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, HighFive::File::ReadOnly};
    validateElectromagDump(sim, hi5);
}

template<typename Simulator, typename Hi5Diagnostic>
void validateParticleDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();

    auto checkParticles = [&](auto& file, auto& particles, auto path) {
        if (!particles.size())
            return;
        std::vector<float> weightV, chargeV, vV;
        file.getDataSet(path + "weight").read(weightV);
        file.getDataSet(path + "charge").read(chargeV);
        file.getDataSet(path + "v").read(vV);
        std::vector<int> iCellV;
        file.getDataSet(path + "iCell").read(iCellV);
        std::vector<float> deltaV;
        file.getDataSet(path + "delta").read(deltaV);

        core::ParticlePacker packer{particles};

        auto first       = packer.empty();
        size_t iCellSize = std::get<2>(first).size();
        size_t deltaSize = std::get<3>(first).size();
        size_t vSize     = std::get<4>(first).size();
        size_t part_idx  = 0;
        while (packer.hasNext())
        {
            auto next = packer.next();

            for (size_t i = 0; i < iCellSize; i++)
                EXPECT_EQ(iCellV[(part_idx * iCellSize) + i], std::get<2>(next)[i]);

            for (size_t i = 0; i < deltaSize; i++)
                EXPECT_FLOAT_EQ(deltaV[(part_idx * deltaSize) + i], std::get<3>(next)[i]);

            for (size_t i = 0; i < vSize; i++)
                EXPECT_FLOAT_EQ(vV[(part_idx * vSize) + i], std::get<4>(next)[i]);

            part_idx++;
        }
    };

    auto checkFile = [&](auto& path, auto tree, auto& particles) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree));
        checkParticles(hifile->file(), particles, path + "/");
    };

    auto visit = [&](GridLayout&, std::string patchID, size_t iLevel) {
        auto path = hi5.getPatchPath(iLevel, patchID);
        for (auto& pop : hybridModel.state.ions)
        {
            checkFile(path, "/ions/pop/" + pop.name() + "/domain", pop.domainParticles());
            checkFile(path, "/ions/pop/" + pop.name() + "/levelGhost", pop.levelGhostParticles());
        }
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


TYPED_TEST(SimulatorTest, particles)
{
    TypeParam sim;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;

    { // scoped to destruct after dump
        Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.particles("/ions/pop/alpha/domain"))
            .addDiagDict(hi5.particles("/ions/pop/alpha/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/alpha/patchGhost"))
            .addDiagDict(hi5.particles("/ions/pop/protons/domain"))
            .addDiagDict(hi5.particles("/ions/pop/protons/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/protons/patchGhost"));
        sim.dump(hi5.dMan);
    }

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, HighFive::File::ReadOnly};
    validateParticleDump(sim, hi5);
}

template<typename Simulator, typename Hi5Diagnostic>
void validateAttributes(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout         = typename Simulator::PHARETypes::GridLayout_t;
    constexpr auto dimension = Simulator::dimension;

    auto& hybridModel = *sim.getHybridModel();
    auto hifile       = hi5.writer.makeFile(hi5.writer.fileString("/EM_B"));

    auto _check_equal = [](auto& group, auto expected, auto key) {
        std::vector<typename decltype(expected)::value_type> attr;
        group.getAttribute(key).read(attr);
        EXPECT_EQ(expected, attr);
    };

    auto visit = [&](GridLayout& grid, std::string patchID, size_t iLevel) {
        auto group = hifile->file().getGroup(hi5.getPatchPath(iLevel, patchID));

        _check_equal(group, grid.origin().toVector(), "origin");
        _check_equal(group, core::Point<uint32_t, dimension>{grid.nbrCells()}.toVector(),
                     "nbrCells");
        _check_equal(group, grid.AMRBox().lower.toVector(), "lower");
        _check_equal(group, grid.AMRBox().upper.toVector(), "upper");
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


TYPED_TEST(SimulatorTest, allFromPython)
{
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    TypeParam sim;
    sim.dump(*sim.dMan);
    sim.dMan.reset();

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.hierarchy;

    Hi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, HighFive::File::ReadOnly};

    validateFluidDump(sim, hi5);
    validateElectromagDump(sim, hi5);
    validateParticleDump(sim, hi5);
    validateAttributes(sim, hi5);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE::SamraiLifeCycle samsam(argc, argv);
    auto ret = RUN_ALL_TESTS();
    return ret;
}
