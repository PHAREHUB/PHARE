
// input/input_1d_ratio_2.txt is unused but a reference

#include <mpi.h>

#include "include.h"

using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;


template<typename Hierarchy, typename HybridModel>
struct FluidHi5Diagnostic : public Hi5Diagnostic<Hierarchy, HybridModel>
{
    FluidHi5Diagnostic(Hierarchy& hierarchy, HybridModel& model, unsigned flags)
        : Hi5Diagnostic<Hierarchy, HybridModel>{hierarchy, model, "fluid", flags}
    {
    }
};

template<typename Simulator, typename Writer, typename ModelView>
void validate_fluid_dump(Simulator& sim, Writer& writer, ModelView& modelView)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();

    auto checkIons = [&](auto& ions, auto patchPath) {
        std::string path(patchPath + "/ions/");

        for (auto& pop : modelView.getIons())
        {
            std::string popPath(path + "pop/" + pop.name() + "/");
            checkField(writer.file(), pop.density(), popPath + "density");
            checkVecField(writer.file(), pop.flux(), popPath + "flux");
        }

        checkField(writer.file(), ions.density(), path + "density");
        checkVecField(writer.file(), ions.velocity(), path + "bulkVelocity");
    };

    auto visit = [&]([[maybe_unused]] GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        auto patchPath = writer.getPatchPath("time", iLevel, patchID);
        checkIons(hybridModel.state.ions, patchPath);
    };

    PHARE::amr::visitHierarchy<GridLayout>(hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.getNumberOfLevels(), hybridModel);
}

TYPED_TEST(SimulatorTest, fluid)
{
    TypeParam sim;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();

    { // Scoped to destruct after dump
        FluidHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.fluid("/ions/density"))
            .addDiagDict(hi5.fluid("/ions/bulkVelocity"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_alpha/density"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_alpha/flux"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_protons/density"))
            .addDiagDict(hi5.fluid("/ions/pop/ions_protons/flux"))
            .dump();
    }

    FluidHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel,
                                                   HighFive::File::ReadOnly};
    validate_fluid_dump(sim, hi5.writer, hi5.modelView);
}


template<typename Hierarchy, typename HybridModel>
struct ElectromagHi5Diagnostic : public Hi5Diagnostic<Hierarchy, HybridModel>
{
    ElectromagHi5Diagnostic(Hierarchy& hierarchy, HybridModel& model, unsigned flags)
        : Hi5Diagnostic<Hierarchy, HybridModel>{hierarchy, model, "electromag", flags}
    {
    }
};


template<typename Simulator, typename Writer>
void validate_electromag_dump(Simulator& sim, Writer& writer)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();

    auto visit = [&]([[maybe_unused]] GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        auto patchPath = writer.getPatchPath("time", iLevel, patchID) + "/";
        auto& B        = hybridModel.state.electromag.B;
        auto& E        = hybridModel.state.electromag.E;
        checkVecField(writer.file(), B, patchPath + B.name());
        checkVecField(writer.file(), E, patchPath + E.name());
    };

    PHARE::amr::visitHierarchy<GridLayout>(hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.getNumberOfLevels(), hybridModel);
}

TYPED_TEST(SimulatorTest, electromag)
{
    TypeParam sim;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();
    { // scoped to destruct after dump
        ElectromagHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.electromag("/EM_B")).addDiagDict(hi5.electromag("/EM_E")).dump();
    }

    ElectromagHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel,
                                                        HighFive::File::ReadOnly};
    validate_electromag_dump(sim, hi5.writer);
}


template<typename Hierarchy, typename HybridModel>
struct ParticlesHi5Diagnostic : public Hi5Diagnostic<Hierarchy, HybridModel>
{
    ParticlesHi5Diagnostic(Hierarchy& hierarchy, HybridModel& model, unsigned flags)
        : Hi5Diagnostic<Hierarchy, HybridModel>{hierarchy, model, "particle", flags}
    {
    }
};

template<typename Simulator, typename Writer>
void validate_particle_dump(Simulator& sim, Writer& writer)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();

    auto checkParticles = [&](auto& particles, auto path) {
        if (!particles.size())
            return;
        std::vector<float> weightV, chargeV, vV;
        writer.file().getDataSet(path + "weight").read(weightV);
        writer.file().getDataSet(path + "charge").read(chargeV);
        writer.file().getDataSet(path + "v").read(vV);
        std::vector<int> iCellV;
        writer.file().getDataSet(path + "iCell").read(iCellV);
        std::vector<float> deltaV;
        writer.file().getDataSet(path + "delta").read(deltaV);

        ParticlePacker packer{particles};

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

    auto visit = [&]([[maybe_unused]] GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        auto patchPath = writer.getPatchPath("time", iLevel, patchID) + "/";
        for (auto& pop : hybridModel.state.ions)
        {
            std::string particlePath(patchPath + "/ions/pop/" + pop.name() + "/");
            checkParticles(pop.domainParticles(), particlePath + "domain/");
            checkParticles(pop.levelGhostParticles(), particlePath + "levelGhost/");
        }
    };

    PHARE::amr::visitHierarchy<GridLayout>(hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.getNumberOfLevels(), hybridModel);
}


TYPED_TEST(SimulatorTest, particles)
{
    TypeParam sim;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();

    { // scoped to destruct after dump
        ParticlesHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel, NEW_HI5_FILE};
        hi5.dMan.addDiagDict(hi5.particles("/ions/pop/ions_alpha/domain"))
            .addDiagDict(hi5.particles("/ions/pop/ions_alpha/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/ions_alpha/patchGhost"))
            .addDiagDict(hi5.particles("/ions/pop/ions_protons/domain"))
            .addDiagDict(hi5.particles("/ions/pop/ions_protons/levelGhost"))
            .addDiagDict(hi5.particles("/ions/pop/ions_protons/patchGhost"))
            .dump();
    }

    ParticlesHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel,
                                                       HighFive::File::ReadOnly};
    validate_particle_dump(sim, hi5.writer);
}

TYPED_TEST(SimulatorTest, attributesRead)
{
    using GridLayout  = typename TypeParam::PHARETypes::GridLayout_t;
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    TypeParam sim;

    auto& hybridModel = *sim.getHybridModel();
    auto& hierarchy   = *sim.getPrivateHierarchy();

    ParticlesHi5Diagnostic<Hierarchy, HybridModel> hi5{hierarchy, hybridModel,
                                                       HighFive::File::ReadOnly};
    auto& writer = hi5.writer;

    auto visit = [&](GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        std::string patchPath = writer.getPatchPath("time", iLevel, patchID) + "/", origin;
        writer.file().getGroup(patchPath).getAttribute("origin").read(origin);
        EXPECT_EQ(gridLayout.origin().str(), origin);
    };

    PHARE::amr::visitHierarchy<GridLayout>(hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.getNumberOfLevels(), hybridModel);
}


TYPED_TEST(SimulatorTest, allFromPython)
{
    using HybridModel = typename TypeParam::HybridModel;
    using Hierarchy   = typename TypeParam::Hierarchy;

    using DiagnosticModelView = typename Hi5Diagnostic<Hierarchy, HybridModel>::DiagnosticModelView;
    using DiagnosticWriter    = typename Hi5Diagnostic<Hierarchy, HybridModel>::DiagnosticWriter;

    TypeParam sim;
    sim.dMan->dump();
    sim.dMan.reset();

    DiagnosticModelView modelView{*sim.getPrivateHierarchy(), *sim.getHybridModel()};
    DiagnosticWriter writer{modelView, "lol.5", HighFive::File::ReadOnly};
    DiagnosticsManager<DiagnosticWriter> dMan{writer};

    validate_fluid_dump(sim, writer, modelView);
    validate_electromag_dump(sim, writer);
    validate_particle_dump(sim, writer);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    PHARE_test::SamraiLifeCycle samsam(argc, argv);
    auto ret = RUN_ALL_TESTS();
    return ret;
}
